#! /usr/bin/env python3
'''
Check each consensus sequence of an MVC output, and report all that have breadth (non-N) below a given threshold
'''

# imports
from pathlib import Path
from sys import argv

# run tool
if __name__ == "__main__":
    if len(argv) != 2 or '-h' in argv or '--help' in argv or '-help' in argv:
        print("USAGE: %s <MVC_output>" % argv[0]); exit(1)
    out_path = Path(argv[1])
    if not out_path.is_dir():
        raise ValueError(f"Directory not found: {argv[1]}")
    print("reference\tlength\tnum_ambiguous\tnum_unambiguous\tbreadth")
    for p in out_path.glob('*.consensus.fas'):
        with open(p, mode='rt') as f:
            seq = f.readlines()[1].strip()
        seq_len = len(seq)
        num_ambig = seq.count('N')
        num_unambig = seq_len - num_ambig
        if seq_len == 0:
            breadth = 0
        else:
            breadth = num_unambig / seq_len
        print(f"{p.name.replace('.consensus.fas','')}\t{seq_len}\t{num_ambig}\t{num_unambig}\t{breadth}")
