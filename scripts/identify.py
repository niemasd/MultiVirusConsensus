#! /usr/bin/env python3
'''
Attempt to identify a list of reference sequences that likely appear in the sample (each with a score), or blank if nothing exists.
'''

# imports
from pathlib import Path
from statistics import mean, stdev
from sys import stderr
from tqdm import tqdm
import argparse

# constants
DEFAULT_THRESH_STDEV = 3
DEFAULT_THRESH_COVERAGE = 1
DEFAULT_THRESH_COMPLETE = 0.1

# run tool
if __name__ == "__main__":
    # parse user args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mvc_output', type=str, help="MVC Output Folder")
    parser.add_argument('--thresh_stdev', type=float, required=False, default=DEFAULT_THRESH_STDEV, help="Number of Standard Deviations Above Mean to Call a Coverage 'Significant'")
    parser.add_argument('--thresh_coverage', type=float, required=False, default=DEFAULT_THRESH_COVERAGE, help="Minimum Possible 'Significant' Coverage Value")
    parser.add_argument('--thresh_complete', type=float, required=False, default=DEFAULT_THRESH_COMPLETE, help="Minimum Possible 'Significant' Completeness Value")
    parser.add_argument('--quiet', action='store_true', help="Suppress Log Messages")
    args = parser.parse_args()
    if args.thresh_stdev <= 0:
        raise ValueError(f"'thresh_stdev' must be positive: {args.thresh_stdev}")
    if args.thresh_coverage < 0:
        raise ValueError(f"'thresh_coverage' must be non-negative: {args.thresh_coverage}")
    if args.thresh_complete < 0 or args.thresh_complete > 1:
        raise ValueError(f"'thresh_complete' must be in the range [0,1]: {args.thresh_complete}")
    args.mvc_output = Path(args.mvc_output)
    if not args.mvc_output.is_dir():
        raise ValueError(f"Directory not found: {args.mvc_output}")

    # load coverages
    cov = dict() # cov[ref] = coverage of ref = total base count in ref / length of ref
    if not args.quiet:
        print(f"Loading coverages from: {args.mvc_output}/*.poscounts.tsv", file=stderr)
    for p in tqdm(list(args.mvc_output.glob('*.poscounts.tsv')), disable=args.quiet):
        with open(p, mode='rt') as f:
            lines = [l.strip() for l in f.readlines()]
        total_count = sum(int(l.split('\t')[-1]) for i, l in enumerate(lines) if i != 0)
        ref_len = int(lines[-1].split('\t')[0])
        if int(lines[1].split('\t')[0]) == 0:
            ref_len += 1
        cov[p.name.replace('.poscounts.tsv','')] = total_count / ref_len

    # load completenesses
    com = dict() # com[ref] = completeness of ref = number of non-N in consensus / consensus length
    if not args.quiet:
        print(f"Loading completenesses from: {args.mvc_output}/*.consensus.fas", file=stderr)
    for p in tqdm(list(args.mvc_output.glob('*.consensus.fas')), disable=args.quiet):
        with open(p, mode='rt') as f:
            seq = ''.join(l.strip() for l in f.readlines()[1:]).upper()
        com[p.name.replace('.consensus.fas','')] = (len(seq) - seq.count('N')) / len(seq)

    # identify significant references
    cov_mean = mean(cov.values())
    cov_std = stdev(cov.values())
    thresh_coverage = max(args.thresh_coverage, cov_mean + (args.thresh_stdev * cov_std))
    cov_sig = {k:v for k, v in cov.items() if v >= thresh_coverage and com[k] >= args.thresh_complete}
    cov_sig_total = sum(cov_sig.values())
    print('score\treference\tcoverage\tcompleteness')
    if len(cov_sig) == 0:
        print('0\tno matches\t0\t0')
    for score, ref in sorted(((v/cov_sig_total, k) for k, v in cov_sig.items()), reverse=True):
        print(f'{score}\t{ref}\t{cov[ref]}\t{com[ref]}')
