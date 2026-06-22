#! /usr/bin/env python3
'''
Check for potential multi-strain/lineage mixtures in an MVC output
'''

# imports
from pathlib import Path
import argparse

# run tool
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mvc_output', type=str, help="MVC Output Folder")
    parser.add_argument('-d', '--min_depth', type=int, default=10, help="Minimum Depth")
    parser.add_argument('-f', '--min_freq', type=float, default=0.9, help="Minimum Consensus Frequency to *Not* be Called Mixture")
    args = parser.parse_args()
    args.mvc_output = Path(args.mvc_output)
    if not args.mvc_output.is_dir():
        raise ValueError(f"Directory not found: {args.mvc_output}")
    for p in args.mvc_output.glob('*.poscounts.tsv'):
        with open(p, mode='rt') as f:
            for row_num, row in enumerate(f):
                if row_num == 0:
                    continue
                pos, a, c, g, t, gap, total = [int(v) for v in row.split('\t')]
                if (total >= args.min_depth) and (max(a, c, g, t, gap) < (args.min_freq * total)):
                    print(f"Possible mixture: {p.name.replace('.poscounts.tsv','')}"); break
