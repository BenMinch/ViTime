# ViTime: Visualizing Virus Timeseries Data
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
3. -r: This is a folder of trimmed fastq reads. They must be paired-end and the paired files should only differ by "_1" and "_2". They also must be zipped with the ".gz" extention.
4. -t: Threads. The more you use the faster it goes. Usually around 10 is sufficient.
5. -minid: This is the Mapping ANI you want to use for generating coverage. This is a number from 0-100 corresponding to percent identity by which to map at. Recommended is 95.
6. -tax: A taxonomy table. This needs to be a tab separated file with at the very least, a column called "query" which has all of your genome names (without extensions) and a column called "order" with the order level taxonomic classification. (TRICK: you can put whatever level classification you want in this order column and it will work fine, you may just need to change the titles of the graphs generated). This file is conviniently generated from my other program PIGv so if you use both this is easy. 
7. -ref: This is the reference csv that has information about your fastq samples. This must include three columns called "SRR_run", "Date", and "Sample_Name". It can also contain columns of environmental variables and numerical values corresponding to them. 
8. -envs: This is where you can list which environmental variables you want to look it. It needs to be a comma separated list that corresponds to columns you have in the reference datasheet.

## Example usage
While I can't include example fastq reads as they are too large for GitHub, provided are some examples of what the reference and taxonomy sheets should look like. Here is an example of what a command will look like.
`python ViTime.py -i genome_folder -o ViTime_Out -r FastQ_folder -t 10 -minid 95 -tax gvclass.tab -ref reference.csv -envs Temperature,Salinity,pH,Location`

## Example outputs
1. **Datasheets**
* **rpkm_matrix_R.csv**: This is a matrix that has a column for genome and order as well as columns for each day of observation representing the RPKM for that day.
* **environment_data.csv**: This just has the raw environmental data for the given days (not very useful).
* **celeb.csv**: This has information about the celebrity score given to each genome (more on this later).
* **Environmental_correlations.csv**: Two csvs representing correlation coefficients and p-values for all environmental variables and total genomes as well as individual clusters (mentioned later)
2. **Graphs**
* **Genome Heatmap**: This is a heatmap of expression (z-score) over time for each genome in your dataset.
![alt text](https://github.com/BenMinch/ViTime/blob/main/images/Genome_heatmap.png)
* **Order abundance**: This is a stacked barplot showing proportion of order level abundance over time.
![alt text](https://github.com/BenMinch/ViTime/blob/main/images/Order_abundance_over_time.png)
* **Genome Clusters**: Genomes are put into 5 clusters based on expression patterns and these cluster memberships are presented in this figure as well as a representative expression profile for one member in that cluster. 
![alt text](https://github.com/BenMinch/ViTime/blob/main/images/Genome_cluster_membership.png)
* **Celebrity Membership**: Genomes are given celebrity scores based on how many times they "boom and bust" over the course of the experiment. Counts of each celebrity score are shown as well as graphs of a representative from each score.
 ![alt text](https://github.com/BenMinch/ViTime/blob/main/images/Celebrity_membership.png)
* **Environmental Correlation barplot**: A plot of correlation coefficients for total correlation and per-cluster correlation with starts for significant p-values. 
 ![alt text](https://github.com/BenMinch/ViTime/blob/main/images/Environmental_correlations_cluster.png)
* **Environmental Correlation linegraphs**: Linegraphs with linear regression lines are plotted for each environmental variable vs abundance.
 ![alt text](https://github.com/BenMinch/ViTime/blob/main/images/Environmental_correlations_total.png)
* **Categories**: Genomes are put into categories based on abundance patterns. A barplot of the number that fall into each category can be found in the Categorize folder.
# Under the Hood/FAQ

### What if I have multiple observations for the same day/timepoint?
The default setting for ViTime is to average all observations taken from the timepoint/day.

### Celebrity Score
The celebrity score is a metric that I developed in order to categorize the "boom and bust" dynamics that are present in some viral populations. It is basically trying to characterize the amount of times a virus population crashes and rises again dramatically (just like a tiktok celebrity). A higher celebrity score corresponds to a higher number of times this happens. The actual formula to obtain +1 to your celebrity score is that the abundance on one day has to change more than +/- 2* CI of the average abundance from the previous observation. 

### Clustering by expression
Genomes are clustered into 5 categories based on an euclidian distance matrix between RPKM observations throughout the timeseries. This clustering can help to parse out cluster or group-specific patterns that are lost when combining all the data together.

### Categorization of Virus genomes
Categorization is done following the protocol found in the following publication:

1. Dart, E., Fuhrman, J. A., & Ahlgren, N. A. (2023). Diverse Marine T4-like Cyanophage Communities Are Primarily Comprised of Low-Abundance Species Including Species with Distinct Seasonal, Persistent, Occasional, or Sporadic Dynamics. Viruses, 15(2), 581.

### What if I want to change RPKM to another metric such as Raw Count values?
This is definitely possible, you just need to tweak the script a little bit. On line 42 of the ViTime.py script, change the -col parameter to 2 instead of 3. 

# Copywrite
ViTime Copyright (C) 2023 Benjamin Minch

This program is free software: you can redistribute it and/or modify it under the terms of the MIT License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
