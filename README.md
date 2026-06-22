# MultiVirusConsensus (MVC)

MultiVirusConsensus (MVC): Fast consensus genome reconstruction of multiple viruses from a mixed sample

# Installation

MVC is written in Python 3 and has the following dependencies:

* [Minimap2](https://github.com/lh3/minimap2)
* [Samtools](https://github.com/samtools/samtools)
* [ViralConsensus](https://github.com/niemasd/ViralConsensus)

We also provide a Docker image with all dependencies installed: [niemasd/multivirusconsensus](https://hub.docker.com/r/niemasd/multivirusconsensus)

# Usage

A help message demonstrating MVC usage can be viewed using the `-h/--help` argument.

## Example

```bash
./MultiVirusConsensus.py -o output -i example/reads.illumina.30X.fq.gz -r example/NC_001802.fas example/NC_045512.fas example/NC_063383.fas -p example/NC_045512.bed
```

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

# Visualizing Results

You can visualize the results of a MultiVirusConsensus run using our [visualization web application](https://niema.net/MultiVirusConsensus/).

# Citing MultiVirusConsensus (MVC)

If you use MVC in your work, please cite:

> Moshiri N (2026). "MultiVirusConsensus: An accurate and efficient open-source pipeline for identification and consensus sequence generation of multiple viruses from mixed samples." *medRxiv*. [doi:10.64898/2026.03.24.26349218](https://doi.org/10.64898/2026.03.24.26349218)
