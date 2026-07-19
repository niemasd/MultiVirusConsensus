#! /usr/bin/env python3
'''
Calculate a summary TSV from a MultiVirusConsensus output
'''

# imports
from csv import writer
from pathlib import Path
from pysam import AlignmentFile
from statistics import mean, median
from subprocess import run
from sys import argv, stdout

# constants
MIN_BASE_QUAL = 30
MIN_COVERAGE = 10
HEADER = [
    'reference',
    'reads_total',
    'bases_total',
    f'bases_q{MIN_BASE_QUAL}',
    f'bases_q{MIN_BASE_QUAL}_prop',
    'reads_mapped',
    'reads_mapped_prop',
    'reference_length',
    f'positions_cov>={MIN_COVERAGE}',
    f'positions_cov>={MIN_COVERAGE}_prop',
    'positions_cov_mean',
    'positions_cov_median',
    'consensus_path',
    'poscounts_path',
    'inscounts_path',
    'mvc_version',
]

# cached values
BAM_STATS = dict()
REF_NAMES = dict()
REF_LENS = dict()
POS_COV_METRICS = dict()
RUN_INFO = dict()

# calculate multiple statistics from a MVC output BAM in a single pass over the alignments
def calc_bam_stats(bam_path):
    # set up dict
    BAM_STATS['reads_unmapped'] = 0
    BAM_STATS['reads_mapped'] = dict() # keys = references, values = counts
    BAM_STATS['bases_total'] = 0
    BAM_STATS[f'bases_q{MIN_BASE_QUAL}'] = 0

    # calculate stats from BAM
    with AlignmentFile(bam_path, 'rb') as bam:
        for aln in bam:
            # count total and qual-passing bases
            BAM_STATS['bases_total'] += len(aln.query_qualities)
            BAM_STATS[f'bases_q{MIN_BASE_QUAL}'] += sum(1 for score in aln.query_qualities if score >= MIN_BASE_QUAL)

            # handle unmapped reads
            if aln.is_unmapped:
                BAM_STATS['reads_unmapped'] += 1

            # handle only primary alignments
            elif not (aln.is_secondary or aln.is_supplementary):
                if aln.reference_name not in BAM_STATS['reads_mapped']:
                    BAM_STATS['reads_mapped'][aln.reference_name] = 0
                BAM_STATS['reads_mapped'][aln.reference_name] += 1

    # calculate final stats
    BAM_STATS['reads_total'] = BAM_STATS['reads_unmapped'] + sum(BAM_STATS['reads_mapped'].values())
    BAM_STATS['reads_mapped_prop'] = {k : v/BAM_STATS['reads_total'] for k, v in BAM_STATS['reads_mapped'].items()}
    BAM_STATS[f'bases_q{MIN_BASE_QUAL}_prop'] = BAM_STATS[f'bases_q{MIN_BASE_QUAL}'] / BAM_STATS['bases_total']

# parse references FASTA
def parse_refs(refs_path):
    REF_NAMES.clear()
    REF_LENS.clear()
    with open(refs_path, 'rt') as ref_f:
        lines = [l.strip() for l in ref_f.read().strip().splitlines()]
    ID = None; seq_len = None
    for line in lines:
        if line.startswith('>'):
            if ID is not None:
                REF_LENS[ID] = seq_len
            ref_name = line[1:]
            ID = ref_name.split()[0].strip()
            seq_len = 0
            REF_NAMES[ID] = ref_name
        else:
            seq_len += len(line)
    if ID is not None:
        REF_LENS[ID] = seq_len
    return REF_LENS

# calculate position coverage metrics
def calc_pos_cov_metrics(out_path):
    # set things up
    POS_COV_METRICS[f'positions_cov>={MIN_COVERAGE}'] = dict()
    POS_COV_METRICS['positions_cov_mean'] = dict()
    POS_COV_METRICS['positions_cov_median'] = dict()

    # calculate metrics from position count TSV files
    for path in out_path.glob('*.poscounts.tsv'):
        ref = path.name.replace('.poscounts.tsv','').strip()
        with open(path, 'rt') as f:
            coverages = [int(line.split()[-1]) for line_num, line in enumerate(f) if line_num != 0]
        POS_COV_METRICS[f'positions_cov>={MIN_COVERAGE}'][ref] = sum(1 for cov in coverages if cov >= MIN_COVERAGE)
        POS_COV_METRICS['positions_cov_mean'][ref] = mean(coverages)
        POS_COV_METRICS['positions_cov_median'][ref] = median(coverages)

    # calculate final metrics
    if len(REF_LENS) == 0:
        parse_refs(out_path / 'references.fas')
    POS_COV_METRICS[f'positions_cov>={MIN_COVERAGE}_prop'] = {k : v/REF_LENS[k] for k, v in POS_COV_METRICS[f'positions_cov>={MIN_COVERAGE}'].items()}

# calculate general info about the run
def calc_run_info(out_path):
    RUN_INFO.clear()
    with open(out_path / 'MultiVirusConsensus.log', 'rt') as f:
        for l in f:
            if "MultiVirusConsensus (MVC) v" in l:
                RUN_INFO['mvc_version'] = l.replace('=','').strip().split()[-1].replace('v','').strip()

# get the value of a given column for a given reference genome
def get_value(out_path, ref, col):
    if len(BAM_STATS) == 0:
        calc_bam_stats(out_path / 'reads.bam')
    if len(REF_LENS) == 0:
        parse_refs(out_path / 'references.fas')
    if len(POS_COV_METRICS) == 0:
        calc_pos_cov_metrics(out_path)
    if len(RUN_INFO) == 0:
        calc_run_info(out_path)
    ref = ref.strip().upper()
    col = col.strip().lower()
    if col == 'reference':
        return REF_NAMES[ref]
    elif col == 'reference_length':
        return REF_LENS[ref]
    elif col == 'consensus_path':
        return out_path / f'{ref}.consensus.fas'
    elif col == 'poscounts_path':
        return out_path / f'{ref}.poscounts.tsv'
    elif col == 'inscounts_path':
        return out_path / f'{ref}.inscounts.json'
    elif col.startswith('positions'):
        return POS_COV_METRICS[col][ref]
    elif col.startswith('reads_mapped'):
        return BAM_STATS[col][ref]
    elif col in RUN_INFO:
        return RUN_INFO[col]
    elif col in BAM_STATS:
        return BAM_STATS[col]
    else:
        raise ValueError(f"Unknown column name: {col}")

# run tool
if __name__ == "__main__":
    if len(argv) != 2 or '-h' in argv or '--help' in argv or '-help' in argv:
        print("USAGE: %s <MVC_output>" % argv[0]); exit(1)
    out_path = Path(argv[1])
    if not out_path.is_dir():
        raise ValueError(f"Directory not found: {argv[1]}")
    refs = [p.name.replace('.consensus.fas','') for p in out_path.glob('*.consensus.fas')]
    tsv = writer(stdout, delimiter='\t')
    tsv.writerow(HEADER)
    for ref in refs:
        tsv.writerow([get_value(out_path, ref, col) for col in HEADER])
