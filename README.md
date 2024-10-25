# MultiVirusConsensus (MVC)

MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses

# Installation

TODO

# Usage

```
usage: MultiVirusConsensus.py [-h] -i READS [READS ...] -r REFERENCE [REFERENCE ...] -o OUTPUT
                              [--quiet] [--threads THREADS] [--include_multimapped]
                              [--minimap2_path MINIMAP2_PATH] [--minimap2_args MINIMAP2_ARGS]
                              [--samtools_path SAMTOOLS_PATH]
                              [--viral_consensus_path VIRAL_CONSENSUS_PATH] [--viral_consensus_args VIRAL_CONSENSUS_ARGS]

MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses

options:
  -h, --help                                                              show this help message and exit
  -i READS [READS ...], --reads READS [READS ...]                         Input Reads (FASTQ)
  -r REFERENCE [REFERENCE ...], --reference REFERENCE [REFERENCE ...]     Reference Genome (FASTA)
  -o OUTPUT, --output OUTPUT                                              Output Folder
  --quiet                                                                 Suppress Log Output
  --threads THREADS                                                       Number of Threads for Minimap2/Samtools
  --include_multimapped                                                   Include Multimapped Reads in Consensus
  --minimap2_path MINIMAP2_PATH                                           Minimap2 Path (default: minimap2)
  --minimap2_args MINIMAP2_ARGS                                           Minimap2 Arguments (default: -x sr)
  --samtools_path SAMTOOLS_PATH                                           Samtools Path (default: samtools)
  --viral_consensus_path VIRAL_CONSENSUS_PATH                             ViralConsensus Path (default: viral_consensus)
  --viral_consensus_args VIRAL_CONSENSUS_ARGS                             ViralConsensus Arguments (default: )
```

## Example

```bash
./MultiVirusConsensus.py -i example/reads.illumina.30X.fq.gz -r example/NC_001802.fas example/NC_045512.fas example/NC_063383.fas -o output
```

## Minimap2 Arguments

By default, we run Minimap2 using its short read preset (`-x sr`). To use a different present, or to manually specify Minimap2 mapping settings, you can provide the Minimap2 arguments you want to use via: `--minimap2_args`

For example, to use the accurate long reads preset (`-x lr:hq`), you could run the following:

```bash
./MultiVirusConsensus.py --minimap2_args '-x lr:hq' -i hq_long_reads.fq.gz -r example/NC_001802.fas example/NC_045512.fas example/NC_063383.fas -o output
```

In general, please only include arguments from the `Indexing:`, `Mapping:`, and `Alignment:` sections of the Minimap2 usage: **do not include arguments from the `Input/Output:` section of the Minimap2 usage!**

## ViralConsensus Arguments

By default, we run ViralConsensus using its default settings. To manually specify ViralConsenus settings, you can provide the ViralConsensus arguments you want to use via: `--viral_consensus_args`

For example, to change the minimum depth setting to 1 (`-d 1`), you could run the following:

```bash
./MultiVirusConsensus.py --viral_consensus_args '-d 1' -i example/reads.illumina.30X.fq.gz -r example/NC_001802.fas example/NC_045512.fas example/NC_063383.fas -o output
```

In general, please only include arguments related to consensus sequence calling: **do not include arguments related to input/output files!**

# Citing MultiVirusConsensus (MVC)

TODO
