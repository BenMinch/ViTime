# ViTime
A reproducable, simple, and powerful pipeline for analyzing viral time-series data. ViTime takes an input of genomes or marker genes and a folder of raw reads corresponding to different sampling days, and outputs dozens of helpful statistics and graphs to see how your populations are changing over time. ViTime will automatically show taxonomy informed patterns, abundance clusters in your data, and correlations with inputted environmental parameters. It is super easy to use and should make timeseries analysis (something that is usually a pain in the butt due to managing a large library of reads) very easy. 

ViTime was made with virus time series data in mind, but here's a secret, it doesn't just work with viruses. As long as you input a list of genomes and proper taxonomy, this tool could be used for literally anything (prokaryotes, eukaryotes, fungi, fish, etc.). 

# Getting Started

### Dependencies
1. Python 3.8  (although should also work on other python 3 builds)
2. Python packages(Numpy, Pandas, Matplotlib)
3. R 4.3
4. R packages(ggplot2 (3.4.2), tidyverse(2.0), pheatmap(1.0.12), tidyr(1.3), gridExtra(2.3))
5. CoverM (0.4.0) (You should download this in a conda environment)

### Installing ViTime

It is as simple as cloning this repository
`git clone https://github.com/BenMinch/ViTime`

# Usage

## Inputs
1. -i: This is a folder with input genomes or marker genes (only one gene/genome per file). They must be in FASTA format
2. -o: This is the name of the output folder you want to store your results in. It can be whatever you want, the program will create it if it doesn't exist.
3. -r: This is a folder of trimmed fastq reads. They must be paired-end and the paired files should only differ by "_1" and "_2"
