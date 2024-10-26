#! /usr/bin/env python3
'''
MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses
'''

# standard imports
from datetime import datetime
from multiprocessing import cpu_count
from os import makedirs
from os.path import isdir, isfile
from subprocess import run
import argparse
import sys

# useful constants
VERSION = '0.0.1'
global QUIET; QUIET = False
global LOGFILE; LOGFILE = None
DEFAULT_NUM_THREADS = cpu_count()
DEFAULT_MINIMAP2_ARGS = '-x sr'
DEFAULT_VIRALCONSENSUS_ARGS = ''

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to log (prefixed by current time)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    if not QUIET:
        print(tmp, file=sys.stderr, end=end); sys.stderr.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# parse user args
def parse_args():
    # parse args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--reads', required=True, type=str, nargs='+', help="Input Reads (FASTQ)")
    parser.add_argument('-r', '--reference', required=True, type=str, nargs='+', help="Reference Genome (FASTA)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Folder")
    parser.add_argument('-bf', '--biobloom_filter', required=False, type=str, default=None, help="BioBloom Filter (for optional host filtering)")
    parser.add_argument('--quiet', action='store_true', help="Suppress Log Output")
    parser.add_argument('--threads', required=False, type=int, default=DEFAULT_NUM_THREADS, help="Number of Threads for Minimap2/Samtools")
    parser.add_argument('--include_multimapped', action='store_true', help="Include Multimapped Reads in Consensus")
    parser.add_argument('--skip_run', action='store_true', help="Skip Running the Analysis Script")
    parser.add_argument('--biobloomcategorizer_path', required=False, type=str, default='biobloomcategorizer', help="BioBloom Categorizer Path")
    parser.add_argument('--minimap2_path', required=False, type=str, default='minimap2', help="Minimap2 Path")
    parser.add_argument('--minimap2_args', required=False, type=str, default=DEFAULT_MINIMAP2_ARGS, help="Minimap2 Arguments")
    parser.add_argument('--samtools_path', required=False, type=str, default='samtools', help="Samtools Path")
    parser.add_argument('--viral_consensus_path', required=False, type=str, default='viral_consensus', help="ViralConsensus Path")
    parser.add_argument('--viral_consensus_args', required=False, type=str, default=DEFAULT_VIRALCONSENSUS_ARGS, help="ViralConsensus Arguments")
    args = parser.parse_args()

    # check args for validity and return
    check_files = args.reads + args.reference
    if args.biobloom_filter is not None:
        if not args.biobloom_filter.lower().endswith('.bf'):
            raise ValueError("BioBloom Filter should have .bf extension: %s" % args.biobloom_filter)
        check_files.append(args.biobloom_filter)
    for fn in check_files:
        if not isfile(fn):
            raise ValueError("File not found: %s" % fn)
    if isdir(args.output) or isfile(args.output):
        raise ValueError("Output exists: %s" % args.output)
    if args.threads < 1:
        raise ValueError("Number of threads must be positive: %d" % args.threads)
    return args

# load FASTA file
def load_fasta(fn):
    if fn.strip().lower().endswith('.gz'):
        infile = gopen(fn, 'rt')
    else:
        infile = open(fn, 'r')
    seqs = dict(); name = None; seq = ''
    for line in infile:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if name is not None:
                if len(seq) == 0:
                    raise ValueError("Malformed FASTA: %s" % fn)
                seqs[name] = seq
            name = l.split()[0][1:]
            if name in seqs:
                raise ValueError("Duplicate sequence ID (%s): %s" % (name, fn))
            seq = ''
        else:
            seq += ''.join(l.split())
    if name is None or len(seq) == 0:
        raise ValueError("Malformed FASTA: %s" % fn)
    seqs[name] = seq
    return seqs

# load reference genomes
def load_references(fns):
    refs = dict()
    for fn in fns:
        curr = load_fasta(fn)
        for k, v in curr.items():
            if k in refs:
                raise ValueError("Duplicate reference ID: %s" % k)
            refs[k] = v
    return refs

# write FASTA from dict where keys are sequence IDs and values are sequences
def write_references(refs, fn):
    if fn.strip().lower().endswith('.gz'):
        f = gopen(fn, 'wt')
    else:
        f = open(fn, 'w')
    for k, v in refs.items():
        f.write('>'); f.write(k); f.write('\n'); f.write(v); f.write('\n')
    f.close()

# write bash script to run analysis
def write_script(
        reads_fns, refs_fn, refs, script_fn,
        include_multimapped=False, threads=DEFAULT_NUM_THREADS,
        biobloom_filter=None, biobloomcategorizer_path='biobloomcategorizer',
        minimap2_path='minimap2', minimap2_args=DEFAULT_MINIMAP2_ARGS,
        samtools_path='samtools',
        viral_consensus_path='viral_consensus', viral_consensus_args=DEFAULT_VIRALCONSENSUS_ARGS):
    out_path = '/'.join(script_fn.split('/')[:-1])
    f = open(script_fn, 'w')
    f.write("#!/usr/bin/env bash\n")
    f.write("# MultiVirusConsensus (MVC) v%s\n" % VERSION)
    f.write("# MVC Command: %s\n" % ' '.join(sys.argv))
    f.write("'%s' -a -t %d %s '%s' " % (minimap2_path, threads, minimap2_args, refs_fn))
    if biobloom_filter is None: # no host filtering (feed Minimap2 the raw reads)
        f.write(' '.join("'%s'" % fn for fn in reads_fns))
    else:                       # host filtering (call BioBloom to filter, and feed Minimap2 the unfiltered reads)
        f.write("<('%s' -c -n -d -t %d -p '%s/biobloom' -f '%s' %s 2> '%s/biobloom.log')" % (biobloomcategorizer_path, threads, out_path, biobloom_filter, ' '.join("'%s'" % fn for fn in reads_fns), out_path))
    f.write(" 2> '%s/minimap2.log'" % out_path)
    f.write(" | tee >('%s' view -@ %d -o '%s/mapped.bam')" % (samtools_path, threads, out_path))
    if not include_multimapped:
        f.write(" | samtools view -h -F 4 -q 1 | tee")
    for ref_ID, ref_seq in refs.items():
        ref_fn = '%s/reference.%s.fas' % (out_path, ref_ID)
        ref_f = open(ref_fn, 'w'); ref_f.write('>%s\n%s\n' % (ref_ID, ref_seq)); ref_f.close()
        f.write(" >(grep -E '^@.+%s\t|\t%s\t' | '%s' -i - -r '%s' -o '%s/%s.consensus.fas' -op '%s/%s.poscounts.tsv' -oi '%s/%s.inscounts.json' %s)" % (ref_ID, ref_ID, viral_consensus_path, ref_fn, out_path, ref_ID, out_path, ref_ID, out_path, ref_ID, viral_consensus_args))
    f.write(" > /dev/null\n")
    f.close()

# main program execution
def main():
    # set things up
    args = parse_args()
    makedirs(args.output)
    global QUIET; QUIET = args.quiet
    global LOGFILE; LOGFILE = open("%s/MultiVirusConsensus.log" % args.output, 'w')
    print_log("=== MultiVirusConsensus (MVC) v%s ===" % VERSION)
    print_log("Command: %s" % ' '.join(sys.argv))
    print_log("Input Viral Reads: %s" % '\t'.join(args.reads))
    print_log("Viral Reference FASTA: %s" % '\t'.join(args.reference))
    print_log("BioBloom Filter: %s" % args.biobloom_filter)
    print_log("Output Directory: %s" % args.output)
    print_log("Include multimapped reads in consensus? %s" % args.include_multimapped)
    print_log("Number of Threads: %s" % args.threads)
    print_log("BioBloom Categorizer Path: %s" % args.biobloomcategorizer_path)
    print_log("Minimap2 Path: %s" % args.minimap2_path)
    print_log("Minimap2 Arguments: %s" % args.minimap2_args)
    print_log("Samtools Path: %s" % args.samtools_path)
    print_log("ViralConsensus Path: %s" % args.viral_consensus_path)
    print_log("ViralConsensus Arguments: %s" % args.viral_consensus_args)

    # merge references into single FASTA
    refs = load_references(args.reference)
    out_refs_path = '%s/references.fas' % args.output
    print_log("Writing merged references FASTA: %s" % out_refs_path)
    write_references(refs, out_refs_path)

    # create and run bash script
    out_script_path = '%s/run.sh' % args.output
    print_log("Writing bash script: %s" % out_script_path)
    write_script(
        args.reads, out_refs_path, refs, out_script_path,
        include_multimapped=args.include_multimapped, threads=args.threads,
        biobloom_filter=args.biobloom_filter, biobloomcategorizer_path=args.biobloomcategorizer_path,
        minimap2_path=args.minimap2_path, minimap2_args=args.minimap2_args,
        samtools_path=args.samtools_path,
        viral_consensus_path=args.viral_consensus_path, viral_consensus_args=args.viral_consensus_args,
    )
    if not args.skip_run:
        print_log("Running bash script...")
        run(['bash', out_script_path])
    print_log("MVC execution complete")
    LOGFILE.close()

# run tool
if __name__ == "__main__":
    main()
