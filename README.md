# MultiVirusConsensus (MVC)

MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses from a mixed sample

# Installation

MVC is written in Python 3 and has the following dependencies:

* `bash`
* [Minimap2](https://github.com/lh3/minimap2)
* [Samtools](https://github.com/samtools/samtools)
* [ViralConsensus](https://github.com/niemasd/ViralConsensus)
* [BioBloom Tools](https://github.com/bcgsc/biobloom) (optional: only used for host-filtering)

We also provide a Docker image with all dependencies installed: [niemasd/multivirusconsensus](https://hub.docker.com/r/niemasd/multivirusconsensus)

# Usage

MVC can be used as follows:

```
usage: MultiVirusConsensus.py [-h] -i READS [READS ...] -r REFERENCE [REFERENCE ...] -o OUTPUT
                              [--quiet] [--threads THREADS] [--include_multimapped]
                              [--minimap2_path MINIMAP2_PATH] [--minimap2_args MINIMAP2_ARGS]
                              [--samtools_path SAMTOOLS_PATH]
                              [--viral_consensus_path VIRAL_CONSENSUS_PATH] [--viral_consensus_args VIRAL_CONSENSUS_ARGS]

MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses

options:
  -h, --help                                                              show this help message and exit
  -i READS [READS ...], --reads READS [READS ...]                         Input Viral Reads (FASTQ)
  -r REFERENCE [REFERENCE ...], --reference REFERENCE [REFERENCE ...]     Reference Viral Genome(s) (FASTA)
  -o OUTPUT, --output OUTPUT                                              Output Folder
  -bf BIOBLOOM_FILTER, --biobloom_filter BIOBLOOM_FILTER                  BioBloom Filter (for optional host filtering)
  --quiet                                                                 Suppress Log Output
  --threads THREADS                                                       Number of Threads for Minimap2/Samtools/BioBloom
  --include_multimapped                                                   Include Multimapped Reads in Consensus
  --skip_run                                                              Skip Running the Analysis Script
  --biobloomcategorizer_path BIOBLOOMCATEGORIZER_PATH                     BioBloom Categorizer Path
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

## Host Filtering (optional)

By default, MVC does not perform host filtering: it assumes that the reads are ready to be processed as-is. However, if you have a sample that needs to be host-filtered, you can provide MVC a [BioBloom](https://github.com/bcgsc/biobloom) filter constructed from the host genome sequence via the optional `-bf/--biobloom_filter` argument. For example, the following command could be used to build a BioBloom filter from the [hg38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) reference human genome using 8 threads:

```bash
biobloommaker -t 8 -p hg38 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

Note that BioBloom supports creating a single filter from a collection of host genome sequences (e.g. multiple different species, or multiple reference sequences from the same species. For more details about how to use BioBloom, as well as how to tweak filter size (and thus memory) vs. accuracy, please refer to the [BioBloom documentation](https://github.com/bcgsc/biobloom?tab=readme-ov-file#2).

## Optional Arguments

While we anticipate that MVC should perform reasonably well using the default settings, some settings can be optionally tweaked as desired.

### Minimap2 Arguments

By default, we run Minimap2 using its short read preset (`-x sr`). To use a different present, or to manually specify Minimap2 mapping settings, you can provide the Minimap2 arguments you want to use via: `--minimap2_args`

For example, to use the accurate long reads preset (`-x lr:hq`), you could run the following:

```bash
./MultiVirusConsensus.py --minimap2_args '-x lr:hq' -i hq_long_reads.fq.gz -r example/NC_001802.fas example/NC_045512.fas example/NC_063383.fas -o output
```

In general, please only include arguments from the `Indexing:`, `Mapping:`, and `Alignment:` sections of the Minimap2 usage: **do not include arguments from the `Input/Output:` section of the Minimap2 usage!**

### ViralConsensus Arguments

By default, we run ViralConsensus using its default settings. To manually specify ViralConsenus settings, you can provide the ViralConsensus arguments you want to use via: `--viral_consensus_args`

For example, to change the minimum depth setting to 1 (`-d 1`), you could run the following:

```bash
./MultiVirusConsensus.py --viral_consensus_args '-d 1' -i example/reads.illumina.30X.fq.gz -r example/NC_001802.fas example/NC_045512.fas example/NC_063383.fas -o output
```

In general, please only include arguments related to consensus sequence calling: **do not include arguments related to input/output files!**

# Citing MultiVirusConsensus (MVC)

If you use MVC in your work, please cite:

> Moshiri N (2024). "MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses from a mixed sample." [https://github.com/niemasd/MultiVirusConsensus](https://github.com/niemasd/MultiVirusConsensus)
