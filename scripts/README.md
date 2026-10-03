This directory contains scripts that may be helpful for post-processing MultiVirusConsensus (MVC) results.

* [`check_mixture.py`](check_mixture.py) attempts to predict if a given sample has a mixture of multiple strains/lineages of any of the included viruses
* [`pairwise_distances.py`](pairwise_distances.py) computes a pairwise distance matrix of consensus sequences in multiple MVC output folders (e.g. to detect cross-contamination)
* [`plot_coverage.py`](plot_coverage.py) produces coverage plots for all reference sequences as a multi-page PDF
* [`summary.py`](summary.py) calculates many summary statistics about the read mapping, base/indel counting, and consensus sequence calling processes and outputs them as a TSV
    * Importantly, it also attempts to predict which viruses likely actually exist in the sample (e.g. for multi-virus surveillance panels)
