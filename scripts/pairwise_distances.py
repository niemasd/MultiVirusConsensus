#! /usr/bin/env python3
'''
Calculate all pairwise distances of all consensus sequences of the same virus from a collection of MultiVirusConsensus output folders.
'''

# imports
from datetime import datetime
from pathlib import Path
from subprocess import run
from tqdm import tqdm
import argparse

# constants
DEFAULT_MIN_COMPLETENESS = 0.1

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mvc_output', nargs='+', type=str, help="MVC Output Folders")
    parser.add_argument('-e', '--email', type=str, required=True, help="Email (for NCBI)")
    parser.add_argument('-o', '--output', type=str, required=True, help="Pairwise Distance Output Folder")
    parser.add_argument('--min_completeness', type=float, required=False, default=DEFAULT_MIN_COMPLETENESS, help="Minimum Completeness to Include Sequence")
    parser.add_argument('--viralmsa_path', type=str, required=False, default='ViralMSA.py', help="Path to 'ViralMSA.py' Executable")
    parser.add_argument('--tn93_path', type=str, required=False, default='tn93', help="Path to 'tn93' Executable")
    args = parser.parse_args()
    args.mvc_output = sorted((Path(v) for v in args.mvc_output), key=lambda x: x.name.strip().lower())
    if len(args.mvc_output) != len({p.name.strip().lower() for p in args.mvc_output}):
        raise ValueError("Duplicate folder names detected. Please ensure unique folder names (not just unique path names)")
    for p in args.mvc_output:
        if not p.is_dir():
            raise ValueError(f"MVC output folder not found: {p}")
        if p.name.strip().lower() == 'reference':
            raise ValueError(f'"reference" is a reserved name; please rename folder: {p}')
    args.output = Path(args.output)
    if args.output.exists():
        raise ValueError(f"Output exists: {args.output}")
    return args

# load consensus sequences: seqs[ref_ID][sample] = consensus sequence
def load_consensus_seqs(mvc_output_paths, min_completeness=DEFAULT_MIN_COMPLETENESS):
    seqs = {ref_ID:dict() for ref_ID in tqdm((p.name.replace('.consensus.fas','') for mvc_output_p in mvc_output_paths for p in mvc_output_p.glob('*.consensus.fas')), desc="Loading Reference IDs")}
    for ref_ID, seq_dict in tqdm(seqs.items(), desc="Loading Consensus Sequences"):
        for mvc_output_p in mvc_output_paths:
            with open(mvc_output_p / f'{ref_ID}.consensus.fas', mode='rt') as f:
                s = ''.join(l.strip() for l in f.read().splitlines()[1:])
                completeness = (len(s) - s.count('N')) / len(s)
                if completeness >= min_completeness:
                    seq_dict[mvc_output_p.name.strip()] = s
    for ref_ID, seq_dict in tqdm(list(seqs.items()), desc="Pruning Empty References"):
        if len(seq_dict) == 0:
            del seqs[ref_ID]
    return seqs

# perform MSA using MAFFT
def run_mafft(seqs, out_dir, mafft_path='mafft'):
    for ref_ID, seq_dict in tqdm(seqs.items(), desc=f"Running: {mafft_path}"):
        fasta_str = f"{'\n'.join(f'>{k}\n{v}' for k, v in seq_dict.items())}\n"
        with open(out_dir / f'{ref_ID}.mafft.aln', mode='wt') as aln_f:
            with open(out_dir / f'{ref_ID}.mafft.log', mode='wt') as log_f:
                run([mafft_path, '--auto', '-'], input=fasta_str, text=True, stdout=aln_f, stderr=log_f, check=True)

# perform MSA using ViralMSA
def run_viralmsa(seqs, email, out_dir, viralmsa_path='ViralMSA.py'):
    for ref_ID, seq_dict in tqdm(seqs.items(), desc=f"Running: {viralmsa_path}"):
        consensus_path = out_dir / f'{ref_ID}.consensus.fas'
        with open(consensus_path, mode='wt') as fas_f:
            for k, v in sorted(seq_dict.items()):
                fas_f.write(f">{k}\n{v}\n")
        run([viralmsa_path, '-q', '-e', email, '-r', ref_ID, '-s', consensus_path, '-o', out_dir / f'{ref_ID}.viralmsa.out'], check=True)

# compute pairwise distances using tn93
def run_tn93(out_dir, tn93_path='tn93'):
    for viralmsa_out_path in tqdm(sorted(out_dir.glob('*.viralmsa.out')), desc=f"Running: {tn93_path}"):
        ref_ID = viralmsa_out_path.name.replace('.viralmsa.out','')
        aln_path = next(viralmsa_out_path.glob('*.aln'))
        with open(out_dir / f'{ref_ID}.tn93.tsv', mode='wt') as tsv_f:
            with open(out_dir / f'{ref_ID}.tn93.log', mode='wt') as log_f:
                run([tn93_path, '-t', '1', '-l', '1', '-a', 'skip', '-D', '\t', aln_path], text=True, stdout=tsv_f, stderr=log_f, check=True)

# load pairwise distances: dists[ref_ID][u][v] = pairwise distance
def load_dists(out_dir):
    dists = dict()
    for tn93_path in tqdm(sorted(out_dir.glob('*.tn93.tsv')), desc="Loading Distances"):
        ref_ID = tn93_path.name.replace('.tn93.tsv','')
        dists[ref_ID] = dict()
        with open(tn93_path, mode='rt') as f:
            for l in f:
                u, v, d = [s.strip() for s in l.split('\t')]
                if not d.endswith('istance'):
                    if u not in dists[ref_ID]:
                        dists[ref_ID][u] = dict()
                    if v not in dists[ref_ID]:
                        dists[ref_ID][v] = dict()
                    dists[ref_ID][u][v] = dists[ref_ID][v][u] = d
    return dists

# write pairwise distances to TSV distance matrix
def write_distance_matrix(dists, out_dir):
    for ref_ID, ref_dists in tqdm(sorted(dists.items()), desc="Writing Distance Matrix"):
        sample_IDs = sorted(ref_dists.keys())
        with open(out_dir / f'{ref_ID}.distance_matrix.tsv', mode='wt') as f:
            f.write(f"\t{'\t'.join(sample_IDs)}\n")
            for i, u in enumerate(sample_IDs):
                f.write(f"{ref_ID}\t")
                for j, v in enumerate(sample_IDs):
                    if i == j:
                        f.write('0')
                    else:
                        f.write(ref_dists[u][v])
                    if j == (len(sample_IDs) - 1):
                        f.write('\n')
                    else:
                        f.write('\t')

# main script logic
def main():
    args = parse_args()
    args.output.mkdir()
    seqs = load_consensus_seqs(args.mvc_output, min_completeness=args.min_completeness)
    run_viralmsa(seqs, args.email, args.output, viralmsa_path=args.viralmsa_path)
    run_tn93(args.output, tn93_path=args.tn93_path)
    dists = load_dists(args.output)
    write_distance_matrix(dists, args.output)

# run script
if __name__ == "__main__":
    main()
