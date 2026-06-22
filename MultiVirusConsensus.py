#! /usr/bin/env python3
'''
MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses
'''

# standard imports
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime
from multiprocessing import cpu_count
from pathlib import Path
from subprocess import DEVNULL, PIPE, STDOUT, Popen, TimeoutExpired, run
from sys import argv, stderr, stdin, stdout
from typing import BinaryIO, Dict, Iterable, List, Optional, Sequence, Tuple
import argparse
import gzip
import re
import shlex

# useful constants
VERSION = '0.1.0'
QUIET = False
LOGFILE = None
KEEP_MULTIMAPPED_OPTIONS = ['all', 'best', 'none']
DEFAULT_NUM_THREADS = cpu_count()
DEFAULT_KEEP_MULTIMAPPED = 'none'
DEFAULT_MINIMAP2_ARGS = '-x sr'
DEFAULT_VIRALCONSENSUS_ARGS = ''
DEFAULT_BUFSIZE = 1048576 # 1 MB
REFS_PROCESS_WARNING_LIMIT = 100
FASTQ_SUFFIXES = {'.fastq', '.fq'}
SAM_BAM_SUFFIXES = {'.sam', '.bam'}
COMPRESSION_SUFFIXES = {'.gz', '.bgz', '.bz2', '.xz', '.zst'}
SAM_RNAME_RE = re.compile(br'^[^\t\r\n]+\t[^\t\r\n]+\t([^\t\r\n]+)\t')
SAM_FLAG_MAPQ_RE = re.compile(br'^[^\t\r\n]+\t([0-9]+)\t[^\t\r\n]*\t[^\t\r\n]*\t([0-9]+)\t')
SAM_HEADER_SN_RE = re.compile(br'(?:^|\t)SN:([^\t\r\n]+)')
SAFE_FILENAME_RE = re.compile(r'[^A-Za-z0-9._-]+')

# open a file
def open_file(path, mode='rt', buffering=DEFAULT_BUFSIZE):
    if path is None:
        if 'w' in mode.lower():
            return stdout
        return stdin
    if isinstance(path, str):
        path = Path(path)
    if path.suffix.strip().lower() == '.gz':
        return gzip.open(path, mode=mode)
    return open(path, mode=mode, buffering=buffering)

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to log (prefixed by current time)
def print_log(s='', end='\n'):
    tmp = f"[{get_time()}] {s}"
    if not QUIET:
        print(tmp, file=stderr, end=end)
        stderr.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end)
        LOGFILE.flush()

# check minimap2
def check_minimap2(minimap2_path='minimap2'):
    command = [str(minimap2_path), '--version']
    try:
        return run(command, capture_output=True, check=True, text=True).stdout.strip()
    except:
        raise RuntimeError(f"Unable to execute dependency: {minimap2_path}")

# check samtools
def check_samtools(samtools_path='samtools'):
    command = [str(samtools_path), 'version']
    try:
        return run(command, capture_output=True, check=True, text=True).stdout.splitlines()[0].strip().split()[-1].strip()
    except:
        raise RuntimeError(f"Unable to execute dependency: {samtools_path}")

# check viral_consensus
def check_viral_consensus(viral_consensus_path='viral_consensus'):
    command = [str(viral_consensus_path), '--version']
    try:
        return run(command, capture_output=True, check=True, text=True).stdout.strip().split()[-1].replace('v','').strip()
    except:
        raise RuntimeError(f"Unable to execute dependency: {viral_consensus_path}")

# parse command-line fragment
def shell_words(s):
    if s is None or str(s).strip() == '':
        return list()
    return shlex.split(str(s))

# get suffix of file path, ignoring compression suffixes
def noncompression_suffix(path):
    suffixes = [s.lower() for s in Path(path).suffixes]
    while suffixes and suffixes[-1] in COMPRESSION_SUFFIXES:
        suffixes.pop()
    if not suffixes:
        return ''
    return suffixes[-1]

# separate read input files as FASTQ or SAM/BAM
def classify_read_inputs(paths):
    fastq_inputs = list()
    alignment_inputs = list()
    for path in paths:
        suffix = noncompression_suffix(path)
        if suffix in FASTQ_SUFFIXES:
            fastq_inputs.append(path)
        elif suffix in SAM_BAM_SUFFIXES:
            alignment_inputs.append(path)
        else:
            raise ValueError(f"Input reads must be FASTQ/SAM/BAM: {path}")
    if not fastq_inputs and not alignment_inputs:
        raise ValueError("At least one FASTQ, SAM, or BAM input is required")
    return fastq_inputs, alignment_inputs

# parse user args
def parse_args():
    # check for -v/--version
    if ('-v' in argv) or ('--version' in argv):
        print(VERSION); exit(0)

    # parse args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--reads', required=True, type=str, nargs='+', help="Input Viral Reads/Alignments (FASTQ/SAM/BAM)")
    parser.add_argument('-r', '--reference', required=True, type=str, nargs='+', help="Reference Viral Genome(s) (FASTA)")
    parser.add_argument('-p', '--primer', required=False, type=str, nargs='*', default=list(), help="Primers (BED)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Folder")
    parser.add_argument('--quiet', action='store_true', help="Suppress Log Output")
    parser.add_argument('--threads', required=False, type=int, default=DEFAULT_NUM_THREADS, help="Number of Threads for Minimap2/Samtools")
    parser.add_argument('--keep_multimapped', required=False, type=str, default=DEFAULT_KEEP_MULTIMAPPED, help=f"What to Keep for Multimapped Reads (options: {', '.join(KEEP_MULTIMAPPED_OPTIONS)})")
    parser.add_argument('--minimap2_path', required=False, type=str, default='minimap2', help="Path to 'minimap2' Executable (for FASTQ input)")
    parser.add_argument('--minimap2_args', required=False, type=str, default=DEFAULT_MINIMAP2_ARGS, help="Arguments for 'minimap2' (for FASTQ input)")
    parser.add_argument('--samtools_path', required=False, type=str, default='samtools', help="Path to 'samtools' Executable")
    parser.add_argument('--viral_consensus_path', required=False, type=str, default='viral_consensus', help="Path to 'viral_consensus' Executable")
    parser.add_argument('--viral_consensus_args', required=False, type=str, default=DEFAULT_VIRALCONSENSUS_ARGS, help="Arguments for 'viral_consensus'")
    parser.add_argument('-v', '--version', action='store_true', help="Print Version Number")
    args = parser.parse_args()

    # convert file/folder names to Path objects
    args.reads = [Path(fn) for fn in args.reads]
    args.reference = [Path(fn) for fn in args.reference]
    args.primer = [Path(fn) for fn in args.primer]
    args.output = Path(args.output)

    # check args for validity and return
    check_paths = args.reads + args.reference + args.primer
    for path in check_paths:
        if not path.is_file():
            raise ValueError(f"File not found: {path}")
    args.fastq_inputs, args.alignment_inputs = classify_read_inputs(args.reads)
    if args.output.exists():
        raise ValueError(f"Output exists: {args.output}")
    if args.threads < 1:
        raise ValueError(f"Number of threads must be positive: {args.threads}")
    args.keep_multimapped = args.keep_multimapped.strip().lower()
    if args.keep_multimapped not in KEEP_MULTIMAPPED_OPTIONS:
        raise ValueError(f"Invalid 'Keep Multimapped' mode ({args.keep_multimapped}). Options: {', '.join(KEEP_MULTIMAPPED_OPTIONS)}")
    return args

# load FASTA file
def load_fasta(path):
    seqs = dict()
    name = None
    seq_chunks = list()
    with open_file(path, 'rt') as infile:
        for line in infile:
            l = line.strip()
            if len(l) == 0:
                continue
            if l[0] == '>':
                if name is not None:
                    seq = ''.join(seq_chunks)
                    if len(seq) == 0:
                        raise ValueError(f"Malformed FASTA: {path}")
                    seqs[name] = seq
                name = l.split()[0][1:]
                if name in seqs:
                    raise ValueError(f"Duplicate sequence ID ({name}): {path}")
                seq_chunks = list()
            else:
                seq_chunks.append(''.join(l.split()))
    seq = ''.join(seq_chunks)
    if name is None or len(seq) == 0:
        raise ValueError(f"Malformed FASTA: {path}")
    seqs[name] = seq
    return seqs

# load BED file
def load_bed(path):
    primers = dict()
    with open_file(path, 'rt') as infile:
        for line in infile:
            if line.strip() == '' or line.startswith('#'):
                continue
            name = line.split('\t')[0].strip()
            if name == '':
                raise ValueError(f"Malformed BED line with empty reference ID in {path}: {line.rstrip()}")
            if name not in primers:
                primers[name] = ''
            primers[name] += line
    return primers

# load reference genomes from FASTA
def load_references(ref_paths, primer_paths=list()):
    # set things up
    refs = dict()
    primers = dict()

    # load reference FASTAs
    for ref_path in ref_paths:
        curr_refs = load_fasta(ref_path)
        for k, v in curr_refs.items():
            if k in refs:
                raise ValueError(f"Duplicate reference ID in reference FASTAs: {k}")
            refs[k] = v

    # load primers
    for bed_path in primer_paths:
        curr_primers = load_bed(bed_path)
        for k, v in curr_primers.items():
            if k not in refs:
                raise ValueError(f"Reference ID ({k}) found in primer BED ({bed_path}) but not in reference FASTA inputs")
            if k in primers:
                raise ValueError(f"Duplicate reference ID in primer BEDs: {k}")
            primers[k] = v

    # return loaded references and primers
    return refs, primers

# convert a single reference sequence to a FASTA record as bytes
def single_reference_to_fasta_bytes(ref_id, ref_seq):
    return f">{ref_id}\n{ref_seq}\n".encode('utf-8')

# convert all reference sequences to a large FASTA file as bytes
def references_to_fasta_bytes(refs):
    chunks = list()
    for ref_id, ref_seq in refs.items():
        chunks.append(single_reference_to_fasta_bytes(ref_id, ref_seq))
    return b''.join(chunks)

# convert a name to a safe name
def safe_basename(name, fallback):
    safe = SAFE_FILENAME_RE.sub('_', name).strip('._-')
    return safe or fallback

# make safe names for references
def make_safe_ref_names(ref_ids):
    mapping = dict()
    used = set()
    for i, ref_id in enumerate(ref_ids, start=1):
        base = safe_basename(ref_id, f'ref_{i}')
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f'{base}_{suffix}'
            suffix += 1
        mapping[ref_id] = candidate
        used.add(candidate)
    return mapping

# helper class for a ViralConsensus job
@dataclass
class ConsensusJob:
    ref_id: str
    safe_name: str
    proc: Popen
    stdin: BinaryIO
    command: List[str]
    reference_path: Path
    primer_path: Optional[Path]
    consensus_path: Path
    poscounts_path: Path
    inscounts_path: Path
    log_path: Path

# helper class for SAM header info
@dataclass
class SamHeaderInfo:
    extra_sq_lines: Dict[bytes, bytes]
    non_sq_header_lines: List[bytes]

# helper class for SAM record routing stats
@dataclass
class SamRoutingStats:
    source_name: str
    num_lines: int = 0
    num_header_lines: int = 0
    num_alignment_lines: int = 0
    routed_alignment_lines: int = 0
    filtered_alignment_lines: int = 0
    unrecognized_or_nonconsensus_alignment_lines: int = 0

# helper class for overall routing stats
@dataclass
class OverallRoutingStats:
    num_lines: int = 0
    num_header_lines: int = 0
    num_alignment_lines: int = 0
    routed_alignment_lines: int = 0
    filtered_alignment_lines: int = 0
    unrecognized_or_nonconsensus_alignment_lines: int = 0

    # add a SAM record's routing stats
    def add(self, stats):
        self.num_lines += stats.num_lines
        self.num_header_lines += stats.num_header_lines
        self.num_alignment_lines += stats.num_alignment_lines
        self.routed_alignment_lines += stats.routed_alignment_lines
        self.filtered_alignment_lines += stats.filtered_alignment_lines
        self.unrecognized_or_nonconsensus_alignment_lines += stats.unrecognized_or_nonconsensus_alignment_lines

# build Minimap2 command
def build_minimap2_command(args, refs_path, num_refs, fastq_paths):
    cmd = [str(args.minimap2_path), '-a', '-t', str(args.threads)]
    cmd.extend(shell_words(args.minimap2_args))
    if args.keep_multimapped == 'all':
        cmd.extend(['--secondary=yes', '-N', str(num_refs)])
    elif args.keep_multimapped == 'best':
        cmd.append('--secondary=no')
    cmd.append(str(refs_path))
    cmd.extend(str(path) for path in fastq_paths)
    return cmd

# build Samtools command for loading an input SAM/BAM
def build_samtools_input_command(args, alignment_path):
    return [str(args.samtools_path), 'view', '-@', str(args.threads), str(alignment_path)]

# build Samtools command to get header from SAM/BAM
def build_samtools_header_command(args, alignment_path):
    return [str(args.samtools_path), 'view', '-H', str(alignment_path)]

# build Samtools command for writing output BAM
def build_samtools_output_command(args, bam_path):
    return [str(args.samtools_path), 'view', '-@', str(args.threads), '-b', '-o', str(bam_path), '-']

# build ViralConsensus command
def build_viral_consensus_command(args, reference_path, primer_path, consensus_path, poscounts_path, inscounts_path):
    cmd = [str(args.viral_consensus_path), '-i', '-', '-r', str(reference_path)]
    if primer_path is not None:
        cmd.extend(['-p', str(primer_path)])
    cmd.extend(['-o', str(consensus_path), '-op', str(poscounts_path), '-oi', str(inscounts_path)])
    cmd.extend(shell_words(args.viral_consensus_args))
    return cmd

# normalize SAM header line
def normalize_sam_header_line(line):
    return line if line.endswith(b'\n') else line + b'\n'

# parse @SQ fields of SAM
def parse_sq_fields(line):
    sn = None
    ln = None
    for field in line.rstrip(b'\r\n').split(b'\t')[1:]:
        if field.startswith(b'SN:'):
            sn = field[3:]
        elif field.startswith(b'LN:'):
            try:
                ln = int(field[3:])
            except ValueError:
                ln = None
    return sn, ln

# collect header metadata from SAM/BAM input
def collect_alignment_header_info(args, alignment_paths, refs):
    ref_lengths_by_bytes = {ref_id.encode('utf-8'): len(seq) for ref_id, seq in refs.items()}
    extra_sq_lines = dict()
    non_sq_header_lines = list()
    seen_non_sq = set()
    for alignment_path in alignment_paths:
        cmd = build_samtools_header_command(args, alignment_path)
        print_log(f"Reading SAM/BAM header: {' '.join(shlex.quote(x) for x in cmd)}")
        completed = run(cmd, stdout=PIPE, stderr=PIPE)
        if completed.returncode != 0:
            stderr_text = completed.stderr.decode('utf-8', errors='replace').strip()
            raise RuntimeError(f"samtools view -H failed for {alignment_path} with exit status {completed.returncode}{': ' + stderr_text if stderr_text else ''}")
        for raw_line in completed.stdout.splitlines(keepends=True):
            line = normalize_sam_header_line(raw_line)
            if line.startswith(b'@HD'):
                continue # MVC emits a single @HD line for the merged stream
            if line.startswith(b'@SQ'):
                sn, ln = parse_sq_fields(line)
                if sn is None:
                    continue
                expected_len = ref_lengths_by_bytes.get(sn)
                if expected_len is not None:
                    if ln is not None and ln != expected_len:
                        ref_name = sn.decode('utf-8', errors='replace')
                        raise ValueError(f"Reference length mismatch for '{ref_name}' in '{alignment_path}': SAM/BAM header LN={ln}, FASTA length={expected_len}. Pre-mapped SAM/BAM inputs should use the same reference sequences supplied with -r/--reference.")
                    continue # the FASTA-derived @SQ line will be emitted instead.
                if ln is None:
                    ref_name = sn.decode('utf-8', errors='replace')
                    raise ValueError(f"Input '{alignment_path}' has @SQ for non-consensus reference '{ref_name}' without a valid LN field")
                existing = extra_sq_lines.get(sn)
                if existing is not None and existing != line:
                    ref_name = sn.decode('utf-8', errors='replace')
                    raise ValueError(f"Conflicting @SQ header lines for non-consensus reference '{ref_name}' across SAM/BAM inputs")
                extra_sq_lines[sn] = line
                continue
            if line.startswith(b'@RG') or line.startswith(b'@CO'):
                if line not in seen_non_sq:
                    non_sq_header_lines.append(line)
                    seen_non_sq.add(line)
    return SamHeaderInfo(extra_sq_lines=extra_sq_lines, non_sq_header_lines=non_sq_header_lines)

# build combined SAM header
def build_combined_sam_header(refs: Dict[str, str], header_info: SamHeaderInfo) -> List[bytes]:
    lines: List[bytes] = [b"@HD\tVN:1.6\tSO:unknown\n"]
    for ref_id, ref_seq in refs.items():
        lines.append(f"@SQ\tSN:{ref_id}\tLN:len(ref_seq)\n".encode('utf-8'))
    for sn in sorted(header_info.extra_sq_lines):
        lines.append(header_info.extra_sq_lines[sn])
    lines.extend(header_info.non_sq_header_lines)
    command_line = ' '.join(argv).replace('\t', ' ').replace('\r', ' ').replace('\n', ' ')
    lines.append(f"@PG\tID:MultiVirusConsensus\tPN:MultiVirusConsensus\tVN:{VERSION}\tCL:{command_line}\n".encode('utf-8'))
    return lines

# write bytes to output file path
def write_bytes_to_path(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)

# start all ViralConsensus jobs
def start_viral_consensus_jobs(args, refs, primers, out_path):
    safe_names = make_safe_ref_names(refs.keys())
    jobs = dict()
    if len(refs) >= REFS_PROCESS_WARNING_LIMIT:
        print_log(f"WARNING: starting {len(refs)} 'viral_consensus' subprocesses at once. This may require higher OS process/open-file limits.")
    for ref_id, ref_seq in refs.items():
        safe = safe_names[ref_id]
        ref_path = out_path / f'{safe}.reference.fas'
        write_bytes_to_path(ref_path, single_reference_to_fasta_bytes(ref_id, ref_seq))
        primer_path = None
        if ref_id in primers:
            primer_path = out_path / f'{safe}.primers.bed'
            write_bytes_to_path(primer_path, primers[ref_id].encode('utf-8'))
        consensus_path = out_path / f'{safe}.consensus.fas'
        poscounts_path = out_path / f'{safe}.poscounts.tsv'
        inscounts_path = out_path / f'{safe}.inscounts.json'
        log_path = out_path / f'{safe}.viral_consensus.log'
        cmd = build_viral_consensus_command(args, ref_path, primer_path, consensus_path, poscounts_path, inscounts_path)
        log_handle = open(log_path, 'wb')
        try:
            proc = Popen(cmd, stdin=PIPE, stdout=log_handle, stderr=STDOUT, bufsize=DEFAULT_BUFSIZE)
        finally:
            log_handle.close()
        if proc.stdin is None:
            raise RuntimeError(f"viral_consensus stdin pipe was not created for reference: {ref_id}")
        jobs[ref_id] = ConsensusJob(ref_id=ref_id, safe_name=safe, proc=proc, stdin=proc.stdin, command=cmd, reference_path=ref_path, primer_path=primer_path, consensus_path=consensus_path, poscounts_path=poscounts_path, inscounts_path=inscounts_path, log_path=log_path)
        print_log(f"Started viral_consensus for reference: {ref_id}")
    return jobs

# parse MAPQ flag
def parse_flag_mapq(line):
    m = SAM_FLAG_MAPQ_RE.match(line)
    if m is None:
        return None
    try:
        return int(m.group(1)), int(m.group(2))
    except ValueError:
        return None

# check if individual alignment/header should be sent to viral_consensus
def consensus_filter_passes(line, keep_multimapped):
    if line.startswith(b'@'):
        return True
    if keep_multimapped == 'all':
        return True
    parsed = parse_flag_mapq(line)
    if parsed is None:
        return False
    flag, mapq = parsed
    if keep_multimapped == 'best':
        return (flag & 256) == 0
    if keep_multimapped == 'none':
        return (flag & 4) == 0 and mapq >= 1
    return True

# check destination of SAM based on header line
def sam_destinations(line, ref_id_by_bytes):
    if line.startswith(b'@'):
        if line.startswith(b'@SQ'):
            m = SAM_HEADER_SN_RE.search(line)
            if m is None:
                return 'drop', None
            ref_id = ref_id_by_bytes.get(m.group(1))
            if ref_id is None:
                return 'drop', None
            return 'ref', ref_id
        return 'all', None # @HD, @RG, @PG, @CO, etc. are valid for every per-reference SAM
    m = SAM_RNAME_RE.match(line)
    if m is None:
        return 'drop', None
    rname = m.group(1)
    if rname == b'*':
        return 'drop', None
    ref_id = ref_id_by_bytes.get(rname)
    if ref_id is None:
        return 'drop', None
    return 'ref', ref_id

# write bytes to stream safely
def safe_write(stream, data, label):
    try:
        stream.write(data)
    except BrokenPipeError as exc:
        raise RuntimeError(f"'{label}' closed its input pipe while MVC was still streaming SAM data") from exc

# close stream
def close_stream(stream, label):
    if stream is None:
        return
    try:
        stream.close()
    except BrokenPipeError:
        print_log(f"'{label}' pipe was already closed")

# terminate processes
def terminate_processes(processes):
    live = [(name, proc) for name, proc in processes if proc is not None and proc.poll() is None]
    for name, proc in live:
        try:
            print_log(f"Terminating '{name}' after error")
            proc.terminate()
        except Exception:
            pass
    for name, proc in live:
        try:
            proc.wait(timeout=5)
        except TimeoutExpired:
            try:
                print_log(f"Killing '{name}' after timeout")
                proc.kill()
            except Exception:
                pass

# wait or raise error
def wait_or_raise(name, proc, log_hint):
    rc = proc.wait()
    if rc != 0:
        suffix = ''
        if log_hint is not None:
            suffix = f" See log: {log_hint}"
        raise RuntimeError(f"'{name}' failed with exit status {rc}.{suffix}")

# route single line of SAM to ViralConsensus
def route_line_to_consensus(line, jobs, ref_id_by_bytes, keep_multimapped, stats=None):
    if not consensus_filter_passes(line, keep_multimapped):
        if stats is not None and not line.startswith(b'@'):
            stats.filtered_alignment_lines += 1
        return
    kind, ref_id = sam_destinations(line, ref_id_by_bytes)
    if kind == 'all':
        for job in jobs.values():
            safe_write(job.stdin, line, f'viral_consensus[{job.ref_id}]')
    elif kind == 'ref' and ref_id is not None:
        job = jobs.get(ref_id)
        if job is not None:
            safe_write(job.stdin, line, f'viral_consensus[{ref_id}]')
            if stats is not None and not line.startswith(b'@'):
                stats.routed_alignment_lines += 1
        elif stats is not None and not line.startswith(b'@'):
            stats.unrecognized_or_nonconsensus_alignment_lines += 1
    elif stats is not None and not line.startswith(b'@'):
        stats.unrecognized_or_nonconsensus_alignment_lines += 1

# output combined SAM header
def emit_combined_header(header_lines, samtools_proc, jobs, refs, keep_multimapped):
    if samtools_proc.stdin is None:
        raise RuntimeError("samtools stdin pipe was not created")
    ref_id_by_bytes = {ref_id.encode('utf-8'): ref_id for ref_id in refs}
    for line in header_lines:
        line = normalize_sam_header_line(line)
        safe_write(samtools_proc.stdin, line, 'samtools')
        route_line_to_consensus(line, jobs, ref_id_by_bytes, keep_multimapped)
    print_log(f"Emitted merged SAM header with {len(header_lines)} line(s)")

# route SAM records from source stream
def route_sam_records_from_stream(source_name, source_stream, samtools_proc, jobs, refs, keep_multimapped):
    ref_id_by_bytes = {ref_id.encode('utf-8'): ref_id for ref_id in refs}
    stats = SamRoutingStats(source_name=source_name)
    if samtools_proc.stdin is None:
        raise RuntimeError(f"samtools stdin pipe was not created")
    for line in iter(source_stream.readline, b''):
        line = normalize_sam_header_line(line)
        stats.num_lines += 1
        if line.strip() == b'':
            continue
        if line.startswith(b'@'):
            stats.num_header_lines += 1
            continue
        stats.num_alignment_lines += 1
        safe_write(samtools_proc.stdin, line, 'samtools')
        route_line_to_consensus(line, jobs, ref_id_by_bytes, keep_multimapped, stats)
    print_log(f"Read {stats.num_lines} SAM line(s) from '{source_name}':")
    print_log(f"- {stats.num_header_lines} header")
    print_log(f"- {stats.num_alignment_lines} alignment")
    print_log(f"- {stats.routed_alignment_lines} routed to viral_consensus")
    print_log(f"- {stats.filtered_alignment_lines} filtered")
    print_log(f"- {stats.unrecognized_or_nonconsensus_alignment_lines} ignored for non-consensus/unknown references")
    return stats

# get manifest lines for jobs
def manifest_lines_for_jobs(jobs):
    lines = ["reference_id\tsafe_name\treference_fasta\tprimer_bed\tconsensus\tposcounts\tinscounts\tlog\n"]
    for job in jobs.values():
        primer_name = f'{job.safe_name}.primers.bed' if job.primer_path is not None else ''
        lines.append(f"{job.ref_id}\t{job.safe_name}\t{job.safe_name}.reference.fas\t{primer_name}\t{job.safe_name}.consensus.fas\t{job.safe_name}.poscounts.tsv\t{job.safe_name}.inscounts.json\t{job.safe_name}.viral_consensus.log\n")
    return lines

# write manifest to output TSV
def write_manifest(out_path, jobs):
    manifest_path = out_path / 'viral_consensus.manifest.tsv'
    manifest_path.write_text(''.join(manifest_lines_for_jobs(jobs)))

# validate ViralConsensus outputs
def validate_viral_consensus_outputs(jobs):
    missing = [path for job in jobs.values() for path in [job.consensus_path, job.poscounts_path, job.inscounts_path] if not path.exists()]
    if missing:
        raise RuntimeError(f"'viral_consensus' completed but expected output file(s) were not created: {', '.join(str(p) for p in missing)}")

# start command line process with log file capturing stderr
def start_process_with_log(cmd, log_path, **popen_kwargs):
    log_f = open(log_path, 'wb')
    try:
        proc = Popen(cmd, stderr=log_f, **popen_kwargs)
    finally:
        log_f.close()
    return proc

# stream alignments from input SAM/BAM
def stream_alignment_inputs(args, alignment_inputs, samtools_proc, jobs, refs, processes, out_path):
    overall = OverallRoutingStats()
    for i, alignment_path in enumerate(alignment_inputs, start=1):
        safe_input_name = safe_basename(alignment_path.name, f'alignment_{i}')
        log_path = out_path / f'samtools_input.{i:03d}.{safe_input_name}.log'
        cmd = build_samtools_input_command(args, alignment_path)
        source_name = f'samtools-input[{alignment_path}]'
        print_log(f"Starting '{source_name}': {' '.join(shlex.quote(x) for x in cmd)}")
        proc = start_process_with_log(cmd, log_path, stdin=None, stdout=PIPE, bufsize=DEFAULT_BUFSIZE)
        processes.append((source_name, proc))
        if proc.stdout is None:
            raise RuntimeError(f"'{source_name} stdout pipe was not created")
        stats = route_sam_records_from_stream(source_name, proc.stdout, samtools_proc, jobs, refs, args.keep_multimapped)
        close_stream(proc.stdout, f'{source_name} stdout')
        wait_or_raise(source_name, proc, log_path)
        overall.add(stats)
    return overall

# stream alignments from Minimap2 results
def stream_minimap2_input(args, refs_path, refs, samtools_proc, jobs, processes, out_path):
    overall = OverallRoutingStats()
    minimap2_log_path = out_path / 'minimap2.log'
    if len(args.fastq_inputs) == 0:
        minimap2_log_path.write_text('No FASTQ inputs were provided; minimap2 was not run.\n')
        print_log('No FASTQ inputs detected; skipping minimap2')
        return overall
    minimap2_cmd = build_minimap2_command(args, refs_path, len(refs), args.fastq_inputs)
    print_log(f"Starting minimap2: {' '.join(shlex.quote(x) for x in minimap2_cmd)}")
    minimap2_proc = start_process_with_log(minimap2_cmd, minimap2_log_path, stdin=None, stdout=PIPE, bufsize=DEFAULT_BUFSIZE)
    processes.append(('minimap2', minimap2_proc))
    if minimap2_proc.stdout is None:
        raise RuntimeError('minimap2 stdout pipe was not created')
    stats = route_sam_records_from_stream('minimap2', minimap2_proc.stdout, samtools_proc, jobs, refs, args.keep_multimapped)
    close_stream(minimap2_proc.stdout, 'minimap2 stdout')
    wait_or_raise('minimap2', minimap2_proc, minimap2_log_path)
    overall.add(stats)
    return overall

# run alignment streaming pipeline
def run_streaming_pipeline(args, refs, primers):
    out_path = args.output
    bam_path = out_path / 'reads.bam'
    samtools_log_path = out_path / 'samtools.log'
    refs_path = out_path / 'references.fas'
    print_log(f"FASTQ input file(s): {'\t'.join(str(path) for path in args.fastq_inputs) if args.fastq_inputs else None}")
    print_log(f"SAM/BAM input file(s): {'\t'.join(str(path) for path in args.alignment_inputs) if args.alignment_inputs else None}")
    print_log(f"Writing merged reference FASTA directly to output folder: {refs_path}")
    write_bytes_to_path(refs_path, references_to_fasta_bytes(refs))
    header_info = collect_alignment_header_info(args, args.alignment_inputs, refs) if args.alignment_inputs else SamHeaderInfo(extra_sq_lines={}, non_sq_header_lines=[])
    combined_header = build_combined_sam_header(refs, header_info)
    processes = list()
    samtools_proc = None
    jobs = dict()
    try:
        jobs = start_viral_consensus_jobs(args, refs, primers, out_path)
        processes.extend((f'viral_consensus[{ref_id}]', job.proc) for ref_id, job in jobs.items())
        write_manifest(out_path, jobs)
        samtools_cmd = build_samtools_output_command(args, bam_path)
        print_log(f"Starting output samtools: {' '.join(shlex.quote(x) for x in samtools_cmd)}")
        samtools_proc = start_process_with_log(samtools_cmd, samtools_log_path, stdin=PIPE, stdout=DEVNULL, bufsize=DEFAULT_BUFSIZE)
        processes.append(('samtools-output', samtools_proc))
        emit_combined_header(combined_header, samtools_proc, jobs, refs, args.keep_multimapped)
        overall = OverallRoutingStats()
        overall.add(stream_alignment_inputs(args, args.alignment_inputs, samtools_proc, jobs, refs, processes, out_path))
        overall.add(stream_minimap2_input(args, refs_path, refs, samtools_proc, jobs, processes, out_path))
        print_log("Overall routed SAM summary:")
        print_log(f"- {overall.num_lines} line(s)")
        print_log(f"- {overall.num_header_lines} source header line(s) skipped")
        print_log(f"- {overall.num_alignment_lines} alignment line(s)")
        print_log(f"- {overall.routed_alignment_lines} routed to viral_consensus")
        print_log(f"- {overall.filtered_alignment_lines} filtered")
        print_log(f"- {overall.unrecognized_or_nonconsensus_alignment_lines} ignored for non-consensus/unknown references")
        close_stream(samtools_proc.stdin if samtools_proc is not None else None, 'samtools stdin')
        for job in jobs.values():
            close_stream(job.stdin, f'viral_consensus[{job.ref_id}] stdin')
        wait_or_raise('samtools-output', samtools_proc, samtools_log_path)
        failed_jobs: List[str] = []
        for ref_id, job in jobs.items():
            rc = job.proc.wait()
            if rc != 0:
                failed_jobs.append(f'{ref_id}(exit={rc})')
        if failed_jobs:
            raise RuntimeError(f"'viral_consensus' failed for: {', '.join(failed_jobs)}. See per-reference logs in: {out_path}")
        validate_viral_consensus_outputs(jobs)
    except Exception:
        terminate_processes(processes)
        raise

# main program logic
def main():
    # parse user args first (so -h works)
    args = parse_args()

    # check dependencies first
    minimap2_version = check_minimap2(args.minimap2_path)
    samtools_version = check_samtools(args.samtools_path)
    viral_consensus_version = check_viral_consensus(args.viral_consensus_path)

    # dependencies exist, so run program
    args.output.mkdir(parents=True)
    global QUIET; QUIET = args.quiet
    global LOGFILE; LOGFILE = open_file(args.output / 'MultiVirusConsensus.log', 'wt')
    try:
        print_log(f"=== MultiVirusConsensus (MVC) v{VERSION} ===")
        print_log(f"Command: {' '.join(argv)}")
        print_log(f"Input Viral Reads/Alignments: {'\t'.join(str(path) for path in args.reads)}")
        print_log(f"FASTQ Input(s): {'\t'.join(str(path) for path in args.fastq_inputs) if args.fastq_inputs else None}")
        print_log(f"SAM/BAM Input(s): {'\t'.join(str(path) for path in args.alignment_inputs) if args.alignment_inputs else None}")
        print_log(f"Reference FASTA(s): {'\t'.join(str(path) for path in args.reference)}")
        print_log(f"Primer BED(s): {'\t'.join(str(path) for path in args.primer) if args.primer else None}")
        print_log(f"Output Directory: {args.output}")
        print_log(f"What to Keep for Multimapped Reads: {args.keep_multimapped}")
        print_log(f"Number of Threads for minimap2/samtools: {args.threads}")
        print_log(f"minimap2 Path: {args.minimap2_path}")
        print_log(f"minimap2 Version: {minimap2_version}")
        print_log(f"minimap2 Arguments: {args.minimap2_args}")
        print_log(f"samtools Path: {args.samtools_path}")
        print_log(f"samtools Version: {samtools_version}")
        print_log(f"viral_consensus Path: {args.viral_consensus_path}")
        print_log(f"viral_consensus Version: {viral_consensus_version}")
        print_log(f"viral_consensus Arguments: {args.viral_consensus_args}")
        refs, primers = load_references(args.reference, args.primer)
        print_log(f"Loaded {len(refs)} reference genome(s)")
        print_log(f"Loaded primer records for {len(primers)} reference genome(s)")
        run_streaming_pipeline(args, refs, primers)
        print_log("MVC execution complete")
    finally:
        if LOGFILE is not None:
            LOGFILE.close()

# run program
if __name__ == '__main__':
    main()
