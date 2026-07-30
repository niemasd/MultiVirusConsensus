#! /usr/bin/env python3
'''
Attempt to identify a list of reference sequences that likely appear in the sample (each with a score), or blank if nothing exists.
'''

# imports
from pathlib import Path
from statistics import mean, stdev
import argparse

# constants
DEFAULT_STD_THRESH = 3
DEFAULT_MIN_THRESH = 1

# run tool
if __name__ == "__main__":
    # parse user args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mvc_output', type=str, help="MVC Output Folder")
    parser.add_argument('-s', '--std_thresh', type=float, required=False, default=DEFAULT_STD_THRESH, help="Number of Standard Deviations Above Mean to Call a Coverage 'Significant'")
    parser.add_argument('-m', '--min_thresh', type=float, required=False, default=DEFAULT_MIN_THRESH, help="Minimum Possible 'Significant' Coverage Value")
    args = parser.parse_args()
    if args.std_thresh <= 0:
        raise ValueError(f"'std_thresh' must be positive: {args.std_thresh}")
    if args.min_thresh < 0:
        raise ValueError(f"'min_thresh' must be non-negative: {args.min_thresh}")
    args.mvc_output = Path(args.mvc_output)
    if not args.mvc_output.is_dir():
        raise ValueError(f"Directory not found: {args.mvc_output}")

    # load data
    cov = dict() # cov[ref] = coverage of ref = total base count in ref / length of ref
    for p in args.mvc_output.glob('*.poscounts.tsv'):
        with open(p, mode='rt') as f:
            lines = [l.strip() for l in f.readlines()]
        total_count = sum(int(l.split('\t')[-1]) for i, l in enumerate(lines) if i != 0)
        ref_len = int(lines[-1].split('\t')[0])
        if int(lines[1].split('\t')[0]) == 0:
            ref_len += 1
        cov[p.name.replace('.poscounts.tsv','')] = total_count / ref_len

    # identify significant references
    cov_mean = mean(cov.values())
    cov_std = stdev(cov.values())
    cov_thresh = max(args.min_thresh, cov_mean + (args.std_thresh * cov_std))
    cov_sig = {k:v for k, v in cov.items() if v >= cov_thresh}
    cov_sig_total = sum(cov_sig.values())
    print('score\treference')
    for k, v in sorted(((v/cov_sig_total, k) for k, v in cov_sig.items()), reverse=True):
        print(f'{k}\t{v}')
