#! /usr/bin/env python3
'''
Create coverage plots from a MultiVirusConsensus output folder
'''

# imports
from glob import glob
from matplotlib.backends.backend_pdf import PdfPages
from os import remove
from os.path import isfile
from pathlib import Path
from seaborn import lineplot
from sys import argv
import matplotlib.pyplot as plt

# run tool
if __name__ == "__main__":
    # parse args
    if len(argv) == 1 or '-h' in argv or '--help' in argv or '-help' in argv:
        print("USAGE: %s <MVC_output> [MVC_output2] [...]" % argv[0]); exit(1)
    mvc_paths = [Path(d) for d in argv[1:]]
    for p in mvc_paths:
        assert p.is_dir(), "MVC output folder not found: %s" % p
        assert len(list(p.glob('*.poscounts.tsv'))) != 0, "No 'poscounts.tsv' files in folder: %s" % p

    # create PDFs
    for mvc_path in mvc_paths:
        fn = str(mvc_path) + '.coverage.pdf'
        if isfile(fn):
            remove(fn)
        with PdfPages(fn) as pdf:
            for tsv in sorted(mvc_path.glob('*.poscounts.tsv')):
                with open(tsv) as tsv_f:
                    counts = [int(l.split('\t')[-1]) for l in tsv_f.readlines()[1:]]
                ax = lineplot(x=list(range(1,len(counts)+1)), y=counts)
                plt.title(tsv.name)
                plt.xlabel("Genome Position")
                plt.ylabel("Coverage")
                plt.xlim(xmin=0, xmax=len(counts))
                plt.ylim(ymin=0)
                plt.tight_layout()
                pdf.savefig(plt.gcf())
                plt.cla(); plt.clf(); plt.close('all')
