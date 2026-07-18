#! /usr/bin/env python3
'''
Calculate a summary TSV from a MultiVirusConsensus output
'''

# imports
from csv import writer
from pathlib import Path
from pysam import AlignmentFile
from subprocess import run
from sys import argv, stdout

# constants
HEADER = ['reference', 'reads_total', 'bases_total', 'bases_q30', 'reads_mapped']
BAM_STATS = dict()

# calculate multiple statistics from a MVC output BAM in a single pass over the alignments
def calc_bam_stats(bam_path):
    # set up dict
    BAM_STATS['reads_unmapped'] = 0
    BAM_STATS['reads_mapped'] = dict() # keys = references, values = counts
    BAM_STATS['bases_total'] = 0
    BAM_STATS['bases_q30'] = 0

    # calc stats
    with AlignmentFile(bam_path, 'rb') as bam:
        for aln in bam:
            # count total and Q30 bases
            BAM_STATS['bases_total'] += len(aln.query_qualities)
            BAM_STATS['bases_q30'] += sum(1 for score in aln.query_qualities if score >= 30)

            # handle unmapped reads
            if aln.is_unmapped:
                BAM_STATS['reads_unmapped'] += 1

            # handle only primary alignments
            elif not (aln.is_secondary or aln.is_supplementary):
                if aln.reference_name not in BAM_STATS['reads_mapped']:
                    BAM_STATS['reads_mapped'][aln.reference_name] = 0
                BAM_STATS['reads_mapped'][aln.reference_name] += 1
    BAM_STATS['reads_total'] = BAM_STATS['reads_unmapped'] + sum(BAM_STATS['reads_mapped'].values())

# get the value of a given column for a given reference genome
def get_value(out_path, ref, col):
    if len(BAM_STATS) == 0:
        calc_bam_stats(out_path / 'reads.bam')
    ref = ref.strip().upper()
    col = col.strip().lower()
    if col == 'reference':
        return ref
    elif col == 'reads_mapped':
        return BAM_STATS['reads_mapped'][ref]
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
