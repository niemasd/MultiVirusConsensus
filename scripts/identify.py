#! /usr/bin/env python3
'''
Attempt to identify a list of reference sequences that likely appear in the sample (each with a score), or blank if nothing exists.
'''

# imports
from pathlib import Path
from statistics import mean, stdev
from sys import argv

# constants
STD_THRESH = 3 # how many standard deviations above the mean to call a coverage 'significant'

# run tool
if __name__ == "__main__":
    if len(argv) != 2 or '-h' in argv or '--help' in argv or '-help' in argv:
        print("USAGE: %s <MVC_output>" % argv[0]); exit(1)
    out_path = Path(argv[1])
    if not out_path.is_dir():
        raise ValueError(f"Directory not found: {argv[1]}")
    cov = dict() # cov[ref] = coverage of ref = total base count in ref / length of ref
    for p in out_path.glob('*.poscounts.tsv'):
        with open(p, mode='rt') as f:
            lines = [l.strip() for l in f.readlines()]
        total_count = sum(int(l.split('\t')[-1]) for i, l in enumerate(lines) if i != 0)
        ref_len = int(lines[-1].split('\t')[0])
        if int(lines[1].split('\t')[0]) == 0:
            ref_len += 1
        cov[p.name.replace('.poscounts.tsv','')] = total_count / ref_len
    cov_mean = mean(cov.values())
    cov_std = stdev(cov.values())
    cov_thresh = cov_mean + (3 * cov_std)
    cov_sig = {k:v for k, v in cov.items() if v >= cov_thresh}
    cov_sig_total = sum(cov_sig.values())
    print('score\treference')
    for k, v in sorted(((v/cov_sig_total, k) for k, v in cov_sig.items()), reverse=True):
        print(f'{k}\t{v}')
