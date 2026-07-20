#! /usr/bin/env python3
'''
Calculate a summary TSV from a MultiVirusConsensus output
'''

# imports
#from Bio.Align import PairwiseAligner
from csv import writer
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path
from pysam import AlignmentFile
from sequence_align.pairwise import needleman_wunsch
from statistics import mean, median
from sys import argv, stderr, stdout
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from warnings import catch_warnings, simplefilter
import argparse

# constants
DEFAULT_MIN_BASE_QUAL = 30
DEFAULT_MIN_COVERAGE = 10

# cached values
BAM_STATS = dict()
REF_NAMES = dict()
REF_LENS = dict()
POS_COV_METRICS = dict()
RUN_INFO = dict()
INDELS = dict()

# get current time
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print log message
def print_log(s='', end='\n', file=stderr):
    print(f"[{get_time()}] {s}", file=file, end=end)

# calculate multiple statistics from a MVC output BAM in a single pass over the alignments
def calc_bam_stats(bam_path, min_base_qual=DEFAULT_MIN_BASE_QUAL):
    # set up dict
    BAM_STATS['reads_unmapped'] = 0
    BAM_STATS['reads_mapped'] = dict() # keys = references, values = counts
    BAM_STATS['bases_total'] = 0
    BAM_STATS[f'bases_q{min_base_qual}'] = 0
    for b in 'acgt':
        BAM_STATS[f'base_{b}'] = 0

    # calculate stats from BAM
    print_log(f"Calculating BAM statistics from: {bam_path}")
    with AlignmentFile(bam_path, 'rb') as bam:
        for aln in tqdm(bam):
            # handle unmapped reads
            if aln.is_unmapped:
                BAM_STATS['reads_unmapped'] += 1

            # handle only primary alignments
            elif not (aln.is_secondary or aln.is_supplementary):
                # handle mapped reads count
                if aln.reference_name not in BAM_STATS['reads_mapped']:
                    BAM_STATS['reads_mapped'][aln.reference_name] = 0
                BAM_STATS['reads_mapped'][aln.reference_name] += 1

                # handle base counts
                BAM_STATS[f'bases_q{min_base_qual}'] += sum(1 for score in aln.query_qualities if score >= min_base_qual)
                BAM_STATS['bases_total'] += len(aln.query_sequence)
                for b in aln.query_sequence.lower():
                    BAM_STATS[f'base_{b}'] += 1

    # calculate final stats
    BAM_STATS['reads_total'] = BAM_STATS['reads_unmapped'] + sum(BAM_STATS['reads_mapped'].values())
    BAM_STATS['reads_mapped_prop'] = {k : v/BAM_STATS['reads_total'] for k, v in BAM_STATS['reads_mapped'].items()}
    BAM_STATS[f'bases_q{min_base_qual}_prop'] = BAM_STATS[f'bases_q{min_base_qual}'] / BAM_STATS['bases_total']
    BAM_STATS['reads_gc_content'] = (BAM_STATS['base_c'] + BAM_STATS['base_g']) / BAM_STATS['bases_total']
    for b in 'acgt':
        BAM_STATS[f'base_{b}_prop'] = BAM_STATS[f'base_{b}'] / BAM_STATS['bases_total']

# parse references FASTA
def parse_refs(refs_path):
    REF_NAMES.clear()
    REF_LENS.clear()
    print_log(f"Parsing references FASTA: {refs_path}")
    with tqdm() as pbar:
        with open(refs_path, 'rt') as ref_f:
            while True:
                tmp = ref_f.readline()
                if not tmp:
                    break
                ref_name = tmp[1:].strip()
                ID = ref_name.split()[0]
                REF_NAMES[ID] = ref_name
                tmp = ref_f.readline()
                if not tmp:
                    raise ValueError(f"Invalid references FASTA: {refs_path}")
                REF_LENS[ID] = len(tmp.strip())
                pbar.update(1)
    return REF_LENS

# calculate position coverage metrics
def calc_pos_cov_metrics(out_path, min_coverage=DEFAULT_MIN_COVERAGE):
    # set things up
    POS_COV_METRICS[f'positions_cov>={min_coverage}'] = dict()
    POS_COV_METRICS['positions_cov_mean'] = dict()
    POS_COV_METRICS['positions_cov_median'] = dict()
    INDELS[f'insertions>={min_coverage}'] = dict()
    INDELS[f'deletions>={min_coverage}'] = dict()

    # calculate metrics from position count TSV files
    print_log("Calculating position coverage metrics...")
    for path in tqdm(list(out_path.glob('*.poscounts.tsv'))):
        ref = path.name.replace('.poscounts.tsv','').strip()
        with open(path, 'rt') as f: # rows are [Position, A, C, G, T, -, Total]
            rows = [[int(x) for x in line.split('\t')] for line_num, line in enumerate(f) if line_num != 0]
        coverages = [row[-1] for row in rows]
        POS_COV_METRICS[f'positions_cov>={min_coverage}'][ref] = sum(1 for cov in coverages if cov >= min_coverage)
        POS_COV_METRICS['positions_cov_mean'][ref] = mean(coverages)
        POS_COV_METRICS['positions_cov_median'][ref] = median(coverages)
        INDELS[f'deletions>={min_coverage}'][ref] = sum(1 for row in rows if row[-2] >= min_coverage)

    # calculate metrics from insertion count JSON files
    print_log("Calculating insertion count metrics...")
    for path in tqdm(list(out_path.glob('*.inscounts.json'))):
        ref = path.name.replace('.inscounts.json','').strip()
        INDELS[f'insertions>={min_coverage}'][ref] = 0 # TODO IMPLEMENT

    # calculate final metrics
    if len(REF_LENS) == 0:
        parse_refs(out_path / 'references.fas')
    POS_COV_METRICS[f'positions_cov>={min_coverage}_prop'] = {k : v/REF_LENS[k] for k, v in POS_COV_METRICS[f'positions_cov>={min_coverage}'].items()}

# calculate general info about the run
def calc_run_info(out_path):
    RUN_INFO.clear()
    with open(out_path / 'MultiVirusConsensus.log', 'rt') as f:
        for l in f:
            if "MultiVirusConsensus (MVC) v" in l:
                RUN_INFO['mvc_version'] = l.replace('=','').strip().split()[-1].replace('v','').strip()

# get the value of a given column for a given reference genome
def get_value(out_path, ref, col, min_base_qual=DEFAULT_MIN_BASE_QUAL, min_coverage=DEFAULT_MIN_COVERAGE):
    if len(BAM_STATS) == 0:
        calc_bam_stats(out_path / 'reads.bam', min_base_qual=min_base_qual)
    if len(REF_LENS) == 0:
        parse_refs(out_path / 'references.fas')
    if len(POS_COV_METRICS) == 0:
        calc_pos_cov_metrics(out_path, min_coverage=min_coverage)
    if len(RUN_INFO) == 0:
        calc_run_info(out_path)
    if len(INDELS) == 0:
        calc_indels(out_path)
    ref = ref.strip().upper()
    col = col.strip().lower()
    if col == 'reference':
        return REF_NAMES[ref]
    elif col == 'reference_length':
        return REF_LENS[ref]
    elif col == 'consensus_path':
        return out_path.resolve() / f'{ref}.consensus.fas'
    elif col == 'poscounts_path':
        return out_path.resolve() / f'{ref}.poscounts.tsv'
    elif col == 'inscounts_path':
        return out_path.resolve() / f'{ref}.inscounts.json'
    elif col.startswith('positions'):
        return POS_COV_METRICS[col][ref]
    elif col.startswith('insertions') or col.startswith('deletions'):
        return INDELS[col][ref]
    elif col.startswith('reads_mapped'):
        return BAM_STATS[col][ref]
    elif col in RUN_INFO:
        return RUN_INFO[col]
    elif col in BAM_STATS:
        return BAM_STATS[col]
    else:
        raise ValueError(f"Unknown column name: {col}")

# main program logic
def main():
    # parse user args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mvc_output', type=str, help="MVC Output Folder")
    parser.add_argument('-o', '--output', type=str, default='stdout', help="Summary TSV Output")
    parser.add_argument('--min_base_qual', type=int, default=DEFAULT_MIN_BASE_QUAL, help="Minimum Base Quality")
    parser.add_argument('--min_coverage', type=int, default=DEFAULT_MIN_COVERAGE, help="Minimum Coverage")
    args = parser.parse_args()
    args.mvc_output = Path(args.mvc_output)
    if not args.mvc_output.is_dir():
        raise ValueError(f"Directory not found: {args.mvc_output}")
    print_log(f"MultiVirusConsensus output folder: {args.mvc_output}")
    args.output = args.output.strip().lower()
    print_log(f"Summary TSV output: {args.output}")
    if args.output == 'stdout':
        args.output = stdout
    else:
        args.output = Path(args.output)
        if args.output.exists():
            raise ValueError(f"Summary output exists: {args.output}")
        args.output = open(args.output, 'wt')

    # define header row
    header = [
        'reference',
        'reads_total',
        'bases_total',
        f'bases_q{args.min_base_qual}',
        f'bases_q{args.min_base_qual}_prop',
        'base_a',
        'base_c',
        'base_g',
        'base_t',
        'base_a_prop',
        'base_c_prop',
        'base_g_prop',
        'base_t_prop',
        'reads_gc_content',
        'reads_mapped',
        'reads_mapped_prop',
        'reference_length',
        f'positions_cov>={args.min_coverage}',
        f'positions_cov>={args.min_coverage}_prop',
        'positions_cov_mean',
        'positions_cov_median',
        f'insertions>={args.min_coverage}',
        f'deletions>={args.min_coverage}',
        'mvc_version',
    ]

    # calculate summary info
    print_log("Loading reference genome IDs...")
    refs = [p.name.replace('.consensus.fas','') for p in args.mvc_output.glob('*.consensus.fas')]
    print_log("Setting up output TSV writer...")
    tsv = writer(args.output, delimiter='\t')
    tsv.writerow(header)
    print_log("Writing summary TSV output...")
    for ref in refs:
        tsv.writerow([get_value(args.mvc_output, ref, col, min_base_qual=args.min_base_qual, min_coverage=args.min_coverage) for col in header])
    args.output.close()

# run tool
if __name__ == "__main__":
    main()
