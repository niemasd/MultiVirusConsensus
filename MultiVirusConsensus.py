#! /usr/bin/env python3
'''
MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses
'''

# standard imports
from datetime import datetime
from multiprocessing import cpu_count
from pathlib import Path
from subprocess import run
from sys import argv, stderr, stdin, stdout
import argparse

# useful constants
VERSION = '0.0.3'
global QUIET; QUIET = False
global LOGFILE; LOGFILE = None
KEEP_MULTIMAPPED_OPTIONS = ['all', 'best', 'none']
DEFAULT_NUM_THREADS = cpu_count()
DEFAULT_KEEP_MULTIMAPPED = 'all'
DEFAULT_MINIMAP2_ARGS = '-x sr'
DEFAULT_VIRALCONSENSUS_ARGS = ''
DEFAULT_BUFSIZE = 1048576 # 1 MB

# open a file
def open_file(path, mode='rt', buffering=DEFAULT_BUFSIZE):
    if path is None:
        if 'r' in mode.lower():
            return stdin
        else:
            return stdout
    elif isinstance(path, str):
        path = Path(path)
    if path.suffix.strip().lower() == '.gz':
        return gopen(path, mode=mode)
    else:
        return open(path, mode=mode, buffering=buffering)

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to log (prefixed by current time)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    if not QUIET:
        print(tmp, file=stderr, end=end); stderr.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# parse user args
def parse_args():
    # parse args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--reads', required=True, type=str, nargs='+', help="Input Viral Reads (FASTQ)")
    parser.add_argument('-r', '--reference', required=True, type=str, nargs='+', help="Reference Viral Genome(s) (FASTA)")
    parser.add_argument('-p', '--primer', required=False, type=str, nargs='*', default=list(), help="Primers (BED)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Folder")
    parser.add_argument('-bf', '--biobloom_filter', required=False, type=str, default=None, help="BioBloom Filter (for optional host filtering)")
    parser.add_argument('--quiet', action='store_true', help="Suppress Log Output")
    parser.add_argument('--threads', required=False, type=int, default=DEFAULT_NUM_THREADS, help="Number of Threads for Minimap2/Samtools/BioBloom")
    parser.add_argument('--keep_multimapped', required=False, type=str, default=DEFAULT_KEEP_MULTIMAPPED, help="What to Keep for Multimapped Reads")
    parser.add_argument('--skip_run', action='store_true', help="Skip Running the Analysis Script")
    parser.add_argument('--biobloomcategorizer_path', required=False, type=str, default='biobloomcategorizer', help="BioBloom Categorizer Path")
    parser.add_argument('--minimap2_path', required=False, type=str, default='minimap2', help="Minimap2 Path")
    parser.add_argument('--minimap2_args', required=False, type=str, default=DEFAULT_MINIMAP2_ARGS, help="Minimap2 Arguments")
    parser.add_argument('--samtools_path', required=False, type=str, default='samtools', help="Samtools Path")
    parser.add_argument('--viral_consensus_path', required=False, type=str, default='viral_consensus', help="ViralConsensus Path")
    parser.add_argument('--viral_consensus_args', required=False, type=str, default=DEFAULT_VIRALCONSENSUS_ARGS, help="ViralConsensus Arguments")
    args = parser.parse_args()

    # convert file/folder names to Path objects
    args.reads = [Path(fn) for fn in args.reads]
    args.reference = [Path(fn) for fn in args.reference]
    args.primer = [Path(fn) for fn in args.primer]
    args.output = Path(args.output)
    if args.biobloom_filter is not None:
        args.biobloom_filter = Path(args.biobloom_filter)
    args.minimap2_path = Path(args.minimap2_path)
    args.samtools_path = Path(args.samtools_path)
    args.viral_consensus_path = Path(args.viral_consensus_path)

    # check args for validity and return
    check_paths = args.reads + args.reference + args.primer
    if args.biobloom_filter is not None:
        if args.biobloom_filter.suffix.strip().lower() != '.bf':
            raise ValueError("BioBloom Filter should have .bf extension: %s" % args.biobloom_filter)
        check_paths.append(args.biobloom_filter)
    for path in check_paths:
        if not path.is_file():
            raise ValueError("File not found: %s" % path)
    reference_stems = {path.stem for path in args.reference}
    for path in args.primer:
        if path.stem not in reference_stems:
            raise ValueError("Primer BED does not have corresponding reference FASTA: %s" % path)
    for path in args.reference:
        try:
            reference_stems.remove(path.stem)
        except KeyError:
            raise ValueError("Duplicate reference genome FASTA filename stem: %s" % path.stem)
    if args.output.exists():
        raise ValueError("Output exists: %s" % args.output)
    if args.threads < 1:
        raise ValueError("Number of threads must be positive: %d" % args.threads)
    args.keep_multimapped = args.keep_multimapped.strip().lower()
    if args.keep_multimapped not in KEEP_MULTIMAPPED_OPTIONS:
        raise ValueError("Invalid 'Keep Multimapped' mode (%s). Options: %s" % (args.keep_multimapped, ', '.join(KEEP_MULTIMAPPED_OPTIONS)))
    return args

# load FASTA file
def load_fasta(path):
    seqs = dict(); name = None; seq = ''
    with open_file(path, 'rt') as infile:
        for line in infile:
            l = line.strip()
            if len(l) == 0:
                continue
            if l[0] == '>':
                if name is not None:
                    if len(seq) == 0:
                        raise ValueError("Malformed FASTA: %s" % path)
                    seqs[name] = seq
                name = l.split()[0][1:]
                if name in seqs:
                    raise ValueError("Duplicate sequence ID (%s): %s" % (name, path))
                seq = ''
            else:
                seq += ''.join(l.split())
    if name is None or len(seq) == 0:
        raise ValueError("Malformed FASTA: %s" % path)
    seqs[name] = seq
    return seqs

# load BED file
def load_bed(path):
    primers = dict()
    with open_file(path, 'rt') as infile:
        for line in infile:
            name = line.split('\t')[0]
            if name not in primers:
                primers[name] = ''
            primers[name] += line
    return primers

# load reference genomes from FASTA
def load_references(ref_paths, primer_paths=list()):
    stem_to_primer_path = {primer_path.stem:primer_path for primer_path in primer_paths}
    refs = dict(); primers = dict()
    for ref_path in ref_paths:
        curr_refs = load_fasta(ref_path)
        for k, v in curr_refs.items():
            if k in refs:
                raise ValueError("Duplicate reference ID in reference FASTAs: %s" % k)
            refs[k] = v
        if ref_path.stem in stem_to_primer_path:
            bed_path = stem_to_primer_path[ref_path.stem]
            curr_primers = load_bed(bed_path)
            for k, v in curr_primers.items():
                if k not in curr_refs:
                    raise ValueError("Reference ID (%s) found in primer BED (%s) but not in reference FASTA (%s)" % (k, bed_path, ref_path))
                if k in primers:
                    raise ValueError("Duplicate reference ID in primer BEDs: %s" % k)
                primers[k] = v
    return refs, primers

# write FASTA from dict where keys are sequence IDs and values are sequences
def write_references(refs, path):
    with open_file(path, 'wt') as f:
        for k, v in refs.items():
            f.write('>'); f.write(k); f.write('\n'); f.write(v); f.write('\n')

# write bash script to run analysis
def write_script(
        reads_paths, refs_path, refs, primers, script_path,
        keep_multimapped='all', threads=DEFAULT_NUM_THREADS,
        biobloom_filter=None, biobloomcategorizer_path='biobloomcategorizer',
        minimap2_path='minimap2', minimap2_args=DEFAULT_MINIMAP2_ARGS,
        samtools_path='samtools',
        viral_consensus_path='viral_consensus', viral_consensus_args=DEFAULT_VIRALCONSENSUS_ARGS):
    out_path = script_path.parent
    with open_file(script_path, 'wt') as f:
        # write script header
        f.write("#!/usr/bin/env bash\n")
        f.write("# MultiVirusConsensus (MVC) v%s\n" % VERSION)
        f.write("# MVC Command: %s\n" % ' '.join(argv))

        # write Minimap2 + BioBloom stuff
        f.write("'%s' -a -t %d %s " % (minimap2_path, threads, minimap2_args))
        if keep_multimapped == 'all':
            f.write("--secondary=yes -N %d " % sum(1 for l in open_file(refs_path,'rt') if l.startswith('>')))
        elif keep_multimapped == 'best':
            f.write("--secondary=no ")
        f.write(" '%s' " % refs_path)
        reads_paths_str = ' '.join("'%s'" % path for path in reads_paths)
        if biobloom_filter is None: # no host filtering (feed Minimap2 the raw reads)
            f.write(reads_paths_str)
        else:                       # host filtering (call BioBloom to filter, and feed Minimap2 the unfiltered reads)
            out_biobloom_path = out_path / 'biobloom'
            biobloom_log_path = out_path / 'biobloom.log'
            f.write("<('%s' -c -n -d -t %d -p '%s' -f '%s' %s 2> '%s')" % (biobloomcategorizer_path, threads, out_biobloom_path, biobloom_filter, reads_paths_str, biobloom_log_path))
        minimap2_log_path = out_path / 'minimap2.log'
        f.write(" 2> '%s'" % minimap2_log_path)

        # write Samtools stuff
        f.write(" | tee >('%s' view -@ %d -o '%s/reads.bam')" % (samtools_path, threads, out_path))
        if keep_multimapped == 'none':
            f.write(" | '%s' view -h -F 4 -q 1 | tee" % samtools_path)

        # write ViralConsensus stuff
        for ref_ID, ref_seq in refs.items():
            ref_path = out_path / ('reference.%s.fas' % ref_ID)
            with open_file(ref_path, 'wt') as ref_f:
                ref_f.write('>%s\n%s\n' % (ref_ID, ref_seq))
            f.write(" >(grep -E '^@.+%s\t|\t%s\t' | '%s' -i - -r '%s'" % (ref_ID, ref_ID, viral_consensus_path, ref_path))
            if ref_ID in primers:
                bed_path = out_path / ('primers.%s.bed' % ref_ID)
                with open_file(bed_path, 'wt') as bed_f:
                    bed_f.write(primers[ref_ID])
                f.write(" -p '%s'" % bed_path)
            consensus_path = out_path / ('%s.consensus.fas' % ref_ID)
            poscounts_path = out_path / ('%s.poscounts.tsv' % ref_ID)
            inscounts_path = out_path / ('%s.inscounts.json' % ref_ID)
            viralconsensus_log_path = out_path / ('viral_consensus.%s.log' % ref_ID)
            f.write(" -o '%s' -op '%s' -oi '%s' %s > '%s' 2>&1)" % (consensus_path, poscounts_path, inscounts_path, viral_consensus_args, viralconsensus_log_path))
        f.write(" > /dev/null\n")

        # make sure the bash script waits for all child processes to finish
        f.write("wait\n")

# main program execution
def main():
    # set things up
    args = parse_args()
    args.output.mkdir(parents=True)
    global QUIET; QUIET = args.quiet
    global LOGFILE; LOGFILE = open_file("%s/MultiVirusConsensus.log" % args.output, 'wt')
    print_log("=== MultiVirusConsensus (MVC) v%s ===" % VERSION)
    print_log("Command: %s" % ' '.join(argv))
    print_log("Input Viral Reads: %s" % '\t'.join(str(path) for path in args.reads))
    print_log("Viral Reference FASTA: %s" % '\t'.join(str(path) for path in args.reference))
    print_log("BioBloom Filter: %s" % args.biobloom_filter)
    print_log("Output Directory: %s" % args.output)
    print_log("What to Keep for Multimapped Reads: %s" % args.keep_multimapped)
    print_log("Number of Threads for Minimap2/Samtools/BioBloom: %s" % args.threads)
    print_log("BioBloom Categorizer Path: %s" % args.biobloomcategorizer_path)
    print_log("Minimap2 Path: %s" % args.minimap2_path)
    print_log("Minimap2 Arguments: %s" % args.minimap2_args)
    print_log("Samtools Path: %s" % args.samtools_path)
    print_log("ViralConsensus Path: %s" % args.viral_consensus_path)
    print_log("ViralConsensus Arguments: %s" % args.viral_consensus_args)

    # merge references into single FASTA
    refs, primers = load_references(args.reference, args.primer)
    out_refs_path = '%s/references.fas' % args.output
    print_log("Writing merged references FASTA: %s" % out_refs_path)
    write_references(refs, out_refs_path)

    # create and run bash script
    out_script_path = args.output / 'run.sh'
    print_log("Writing bash script: %s" % out_script_path)
    write_script(
        args.reads, out_refs_path, refs, primers, out_script_path,
        keep_multimapped=args.keep_multimapped, threads=args.threads,
        biobloom_filter=args.biobloom_filter, biobloomcategorizer_path=args.biobloomcategorizer_path,
        minimap2_path=args.minimap2_path, minimap2_args=args.minimap2_args,
        samtools_path=args.samtools_path,
        viral_consensus_path=args.viral_consensus_path, viral_consensus_args=args.viral_consensus_args,
    )
    if not args.skip_run:
        print_log("Running bash script...")
        with open_file(out_script_path.with_suffix('.stdout.log'), 'wt') as run_stdout:
            with open_file(out_script_path.with_suffix('.stderr.log'), 'wt') as run_stderr:
                run(['bash', out_script_path], stdout=run_stdout, stderr=run_stderr)
    print_log("MVC execution complete")
    LOGFILE.close()

# run tool
if __name__ == "__main__":
    main()
