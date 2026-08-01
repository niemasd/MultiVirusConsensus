#! /usr/bin/env python3
'''
Calculate a summary TSV from a MultiVirusConsensus output.

To run on multiple samples `SAMPLE.summary.tsv` and prepend a "sample" column:

(head -1 "$(ls *.summary.tsv | head -1)" | sed 's/^/sample\t/' && for f in *.summary.tsv ; do cat "$f" | grep -v 'reference' | sed "s/^/$(echo $f | cut -d'.' -f1)\t/" ; done ) > summary.tsv
'''

# imports
from datetime import datetime
from json import load as jload
from pandas import DataFrame, read_csv, Series
from pathlib import Path
from pysam import AlignmentFile
from statistics import quantiles
from sys import stderr
from tqdm import tqdm
import argparse

# constants
DEFAULT_FILTER_COVERAGE = 1
DEFAULT_FILTER_COMPLETENESS = 0.1

# get current time
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print log message
def print_log(s='', end='\n', file=stderr, quiet=False):
    if not quiet:
        print(f"[{get_time()}] {s}", file=file, end=end)

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mvc_output', type=str, help="MVC Output Folder")
    parser.add_argument('-o', '--output', type=str, required=False, default='stdout', help="Summary TSV Output")
    parser.add_argument('--filter_coverage', type=float, required=False, default=DEFAULT_FILTER_COVERAGE, help="Filter: Minimum Fold-Coverage")
    parser.add_argument('--filter_completeness', type=float, required=False, default=DEFAULT_FILTER_COMPLETENESS, help="Filter: Minimum Consensus Completeness")
    parser.add_argument('--quiet', action='store_true', help="Suppress Log Messages")
    args = parser.parse_args()
    args.mvc_output = Path(args.mvc_output)
    if not args.mvc_output.is_dir():
        raise ValueError(f"Directory not found: {args.mvc_output}")
    args.output = args.output.strip()
    if args.output != 'stdout':
        args.output = Path(args.output)
        if args.output.exists():
            raise ValueError(f"Output exists: {args.output}")
    return args

# load MultiVirusConsensus log
def load_mvc_log(path, quiet=False):
    print_log(f"Loading MVC log information: {path}", quiet=quiet)
    meta = dict()
    with open(path, mode='rt') as f:
        for line in f:
            if 'MultiVirusConsensus (MVC) v' in line:
                meta['version'] = line.split('MultiVirusConsensus (MVC) v')[1].split()[0]
    return meta

# load reference sequences
def load_ref_seqs(path, suffix='.reference.fas', quiet=False):
    print_log(f"Loading reference sequences from: {path / ('*' + suffix)}", quiet=quiet)
    ref_seqs = dict()
    for p in tqdm(list(path.glob(f'*{suffix}')), disable=quiet):
        ref = p.name.replace(suffix, '')
        with open(p, mode='rt') as f:
            f.readline()
            ref_seqs[ref] = ''.join(l.strip() for l in f.readlines())
    return ref_seqs

# parse ViralConsensus consensus sequence header line
def parse_viral_consensus_header(header):
    header = header.replace('(',' ').replace(')',' ')
    return {
        'version': header.split('viral_consensus v')[1].split()[0].strip(),
        'min_qual': int(header.split('--min_qual ')[1].split()[0].strip()),
        'min_depth': int(header.split('--min_depth ')[1].split()[0].strip()),
        'min_freq': float(header.split('--min_freq ')[1].split()[0].strip()),
    }

# load consensus sequences + ViralConsensus arguments
def load_consensus_seqs(path, suffix='.consensus.fas', quiet=False):
    print_log(f"Loading consensus sequences from: {path / ('*' + suffix)}", quiet=quiet)
    consensus_seqs = dict()
    viral_consensus_meta = None
    for p in tqdm(list(path.glob(f'*{suffix}')), disable=quiet):
        ref = p.name.replace(suffix, '')
        with open(p, mode='rt') as f:
            header = f.readline()
            if viral_consensus_meta is None:
                viral_consensus_meta = parse_viral_consensus_header(header)
            consensus_seqs[ref] = ''.join(l.strip() for l in f.readlines())
    return consensus_seqs, viral_consensus_meta

# load insertion counts
def load_ins_counts(path, suffix='.inscounts.json', quiet=False):
    print_log(f"Loading insertion counts from: {path / ('*' + suffix)}", quiet=quiet)
    ins_counts = dict()
    for p in tqdm(list(path.glob(f'*{suffix}')), disable=quiet):
        ref = p.name.replace(suffix, '')
        with open(p, mode='rt') as f:
            ins_counts[ref] = jload(f)
    return ins_counts

# load position counts
def load_pos_counts(path, suffix='.poscounts.tsv', quiet=False):
    print_log(f"Loading position counts from: {path / ('*' + suffix)}", quiet=quiet)
    return {p.name.replace(suffix,'') : read_csv(p,sep='\t').rename(columns=str.upper) for p in tqdm(list(path.glob(f'*{suffix}')), disable=quiet)}

# calculate statistics from MVC output BAM
def calc_bam_stats(bam_path, min_base_qual, quiet=False):
    # set things up
    print_log(f"Calculating BAM statistics from: {bam_path}", quiet=quiet)
    bam_stats = { # dicts have keys = reference IDs
        'reads_unmapped': 0,
        'reads_mapped': dict(),
        'bases_total': 0,
        f'bases_q{min_base_qual}': 0,
        'base_a': 0,
        'base_c': 0,
        'base_g': 0,
        'base_t': 0,
        'insert_25th': dict(),
        'insert_median': dict(),
        'insert_75th': dict(),
    }
    r1_start_end = dict() # r1_start_end[ref][read] = (start, end)
    r2_start_end = dict() # r2_start_end[ref][read] = (start, end)

    # calculate stats from BAM
    with AlignmentFile(bam_path, 'rb') as bam:
        for aln in tqdm(bam, disable=quiet):
            # handle base counts
            if aln.is_unmapped or not (aln.is_secondary or aln.is_supplementary):
                bam_stats[f'bases_q{min_base_qual}'] += sum(1 for score in aln.query_qualities if score >= min_base_qual)
                bam_stats['bases_total'] += len(aln.query_sequence)
                for b in aln.query_sequence.lower():
                    bam_stats[f'base_{b}'] += 1

            # handle unmapped reads
            if aln.is_unmapped:
                bam_stats['reads_unmapped'] += 1

            # handle only primary alignments
            elif not (aln.is_secondary or aln.is_supplementary):
                # handle mapped reads count
                if aln.reference_name not in bam_stats['reads_mapped']:
                    bam_stats['reads_mapped'][aln.reference_name] = 0
                bam_stats['reads_mapped'][aln.reference_name] += 1

                # keep track of details needed for insert size calculationa
                if aln.reference_name not in r1_start_end:
                    r1_start_end[aln.reference_name] = dict()
                    r2_start_end[aln.reference_name] = dict()
                if aln.is_read1:
                    r1_start_end[aln.reference_name][aln.query_name] = (aln.reference_start, aln.reference_end)
                elif aln.is_read2:
                    r2_start_end[aln.reference_name][aln.query_name] = (aln.reference_start, aln.reference_end)

    # calculate final stats
    bam_stats['reads_total'] = bam_stats['reads_unmapped'] + sum(bam_stats['reads_mapped'].values())
    if bam_stats['reads_total'] == 0:
        bam_stats['reads_mapped_prop'] = {k:0 for k in bam_stats['reads_mapped']}
    else:
        bam_stats['reads_mapped_prop'] = {k : v/bam_stats['reads_total'] for k, v in bam_stats['reads_mapped'].items()}
    if bam_stats['bases_total'] == 0:
        bam_stats[f'bases_q{min_base_qual}_prop'] = 0
        bam_stats['reads_gc_content'] = 0
        for b in 'acgt':
            bam_stats[f'base_{b}_prop'] = 0
    else:
        bam_stats[f'bases_q{min_base_qual}_prop'] = bam_stats[f'bases_q{min_base_qual}'] / bam_stats['bases_total']
        bam_stats['reads_gc_content'] = (bam_stats['base_c'] + bam_stats['base_g']) / bam_stats['bases_total']
        for b in 'acgt':
            bam_stats[f'base_{b}_prop'] = bam_stats[f'base_{b}'] / bam_stats['bases_total']
    for ref in r1_start_end:
        insert_sizes = list()
        for k in r1_start_end[ref]:
            if k in r2_start_end[ref]:
                r1_start, r1_end = r1_start_end[ref][k]
                r2_start, r2_end = r2_start_end[ref][k]
                insert_sizes.append(max(r1_end, r2_end) - min(r1_start, r2_start))
        if len(insert_sizes) == 0:
            bam_stats['insert_25th'][ref], bam_stats['insert_median'][ref], bam_stats['insert_75th'][ref] = 0, 0, 0
        elif len(insert_sizes) == 1:
            bam_stats['insert_25th'][ref], bam_stats['insert_median'][ref], bam_stats['insert_75th'][ref] = insert_sizes[0], insert_sizes[0], insert_sizes[0]
        else:
            bam_stats['insert_25th'][ref], bam_stats['insert_median'][ref], bam_stats['insert_75th'][ref] = quantiles(insert_sizes, n=4)
    return bam_stats

# calculate coverage statistics
def calc_coverage_stats(pos_counts, ins_counts, ref_seqs, min_coverage, quiet=False):
    print_log(f"Calculating coverage statistics from position counts...", quiet=quiet)
    coverage_stats = dict()
    base_cols = ['A', 'C', 'G', 'T']
    allele_cols = base_cols + ['-']
    for ref, df in tqdm(pos_counts.items(), disable=quiet):
        coverage_stats[ref] = dict()
        coverage_stats[ref][f'positions_cov>={min_coverage}'] = (df['TOTAL'] >= min_coverage).sum()
        coverage_stats[ref][f'positions_cov>={min_coverage}_prop'] = coverage_stats[ref][f'positions_cov>={min_coverage}'] / len(df)
        coverage_stats[ref]['positions_cov_mean'] = df['TOTAL'].mean()
        coverage_stats[ref]['positions_cov_median'] = df['TOTAL'].median()
        allele_counts = df[allele_cols]
        max_allele_count = allele_counts.max(axis=1)
        winner = allele_counts.idxmax(axis=1)
        allele_passes_coverage = max_allele_count.ge(min_coverage)
        reference_base = Series(list(ref_seqs[ref].upper()), index=df.index)
        coverage_stats[ref][f'deletions>={min_coverage}'] = int((allele_passes_coverage & winner.eq('-')).sum())
        coverage_stats[ref][f'substitutions>={min_coverage}'] = int((allele_passes_coverage & winner.isin(base_cols) & winner.ne(reference_base)).sum())
    print_log(f"Calculating coverage statistics from insertion counts...", quiet=quiet)
    for ref, count_dicts in tqdm(ins_counts.items(), disable=quiet):
        coverage_stats[ref][f'insertions>={min_coverage}'] = sum(1 for count_dict in count_dicts.values() if (('' not in count_dict) or (count_dict[''] != max(count_dict.values()))) and (max(count_dict.values()) >= min_coverage))
    print_log(f"Calculating matches between consensus and reference sequences...", quiet=quiet)
    for ref, ref_seq in tqdm(ref_seqs.items(), disable=quiet):
        compare_length = coverage_stats[ref][f'positions_cov>={min_coverage}'] + coverage_stats[ref][f'insertions>={min_coverage}']
        if compare_length == 0:
            coverage_stats[ref][f'reference_matches>={min_coverage}'] = 0
            coverage_stats[ref][f'reference_matches>={min_coverage}_prop'] = 0
        else:
            coverage_stats[ref][f'reference_matches>={min_coverage}'] = compare_length - coverage_stats[ref][f'insertions>={min_coverage}'] - coverage_stats[ref][f'deletions>={min_coverage}'] - coverage_stats[ref][f'substitutions>={min_coverage}']
            coverage_stats[ref][f'reference_matches>={min_coverage}_prop'] = coverage_stats[ref][f'reference_matches>={min_coverage}'] / compare_length
    return coverage_stats

# main script logic
def main():
    # load user data
    args = parse_args()
    mvc_meta = load_mvc_log(args.mvc_output / 'MultiVirusConsensus.log', quiet=args.quiet)
    ref_seqs = load_ref_seqs(args.mvc_output, quiet=args.quiet)
    refs = sorted(ref_seqs.keys())
    consensus_seqs, viral_consensus_meta = load_consensus_seqs(args.mvc_output, quiet=args.quiet)
    pos_counts = load_pos_counts(args.mvc_output, quiet=args.quiet)
    ins_counts = load_ins_counts(args.mvc_output, quiet=args.quiet)

    # calculate statistics
    min_base_qual = viral_consensus_meta['min_qual']
    min_coverage = viral_consensus_meta['min_depth']
    bam_stats = calc_bam_stats(args.mvc_output / 'reads.bam', min_base_qual, quiet=args.quiet)
    coverage_stats = calc_coverage_stats(pos_counts, ins_counts, ref_seqs, min_coverage, quiet=args.quiet)

    # fill in missing values for convenience
    for k, v in bam_stats.items():
        if isinstance(v, dict):
            for ref in refs:
                if ref not in v:
                    v[ref] = 0

    # build dataframe
    data = {
        'mvc_version': [mvc_meta['version']]*len(refs),
        'viral_consensus_version': [viral_consensus_meta['version']]*len(refs),
        'reference': refs,
        'reference_length': [len(ref_seqs[ref]) for ref in refs],
    }
    for k in ['reads_total', 'bases_total', f'bases_q{min_base_qual}', f'bases_q{min_base_qual}_prop', 'reads_gc_content']:
        data[k] = [bam_stats[k]]*len(refs)
    for k in ['reads_mapped', 'reads_mapped_prop', 'insert_25th', 'insert_median', 'insert_75th']:
        data[k] = [bam_stats[k][ref] for ref in refs]
    for k in [f'positions_cov>={min_coverage}', f'positions_cov>={min_coverage}_prop', 'positions_cov_mean', 'positions_cov_median', f'substitutions>={min_coverage}', f'insertions>={min_coverage}', f'deletions>={min_coverage}', f'reference_matches>={min_coverage}', f'reference_matches>={min_coverage}_prop']:
        data[k] = [coverage_stats[ref][k] for ref in refs]
    df = DataFrame(data)
    likely_exists_key = f'likely_exists(fold_coverage>={args.filter_coverage};completeness>={args.filter_completeness})'
    df[likely_exists_key] = (df['positions_cov_mean'] >= args.filter_coverage) & (df[f'positions_cov>={min_coverage}_prop'] >= args.filter_completeness)
    df = df.sort_values(by=[likely_exists_key, 'positions_cov_mean', f'positions_cov>={min_coverage}_prop'], ascending=[False, False, False]).reset_index(drop=True)

    # write dataframe to output as TSV
    if args.output == 'stdout':
        from sys import stdout as out_f
    elif args.output.suffix.lower() == '.gz':
        out_f = gopen(args.output, 'wt')
    else:
        out_f = open(args.output, 'wt')
    df.to_csv(out_f, sep='\t', index=False)

# run script
if __name__ == "__main__":
    main()
