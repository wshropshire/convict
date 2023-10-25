# COpy Number Variant quantifICation Tool (CONVICT) 

The purpose of CONVICT is to estimate gene copy number variants of interest using short- or long-reads. This is accomplished through querying a target gene database of interest and calculating the ratio of coverage depths between the aforementioned target gene hits and a user provider set of housekeeping and/or control genes (e.g. cgMLST schema). Regions of high variance in coverage depths of open reading frames (e.g. 5' and 3' ends) are controlled through a binning approach that measures coefficient of variation (CV) per binned ORF region and removes regions iteratively that contribute the greatest variance in coverage depth (CV>0.5) relative to overall mean ORF coverage depth. The default databases currently used are ResFinder (Accessed 2021-11-09) and the KmerFinder_db; users have the option to provide an optional target or species database. KmerFinder_db must be downloaded through the [Center for Genomic Epidemiology (CGE)](https://bitbucket.org/genomicepidemiology/kma/src/master/) ftp server with instructions on their main branch. The primary output is an estimate of normalized coverage depths in sample_convict_results.tsv. 

## Author

[William Shropshire](https://twitter.com/The_Real_Shrops)

## Dependencies

These are the current dependencies that I have installed through conda to run the convict.py script:
* python (3.9.7)
* bowtie2 (2.4.5)
* samtools (1.14)
* minimap2 (2.24) - For Oxford Nanopore Technologies (ONT) long-read alignment
* bbmap (38.93)
* bcftools (1.14)
* bedtools (2.30.0)
and Python3 packages:
* biopython (1.79)
* pandas (1.4.0)
* numpy (1.22.1)

You can test other versions, but there is no guarantee that they will behave well with convict.py 
My future to-do list is create a docker container and/or a bespoke Conda package.

kma (KMA-1.3.24a) and kmerresistance (KmerResistance-2.2.0) can be downloaded with instructions [here](https://bitbucket.org/genomicepidemiology/kma/src/master/) and [here](https://bitbucket.org/genomicepidemiology/kmerresistance/src/master/) respectively. 

This GitHub repository includes ResFinder database that has already been KMA indexed. If using a customary target gene database or species database, you must though use kma index first on your nucleotide fasta file and/or target database of interest: 

```
kma index -i <species_db.fasta> -o <species_db> -Sparse ATG
```

Note that the kma species database needs to be sparse indexed for kmerresistance to properly work. Please indicate the absolute pathway of this indexed database in the convict python script. 

## Usage
An example usage with short reads:
```
$python3 convict.py -t 2 -s sample_name -o outdir -R1 read1.fastq.gz -R2 read2.fastq.gz -c control_genes.fasta
```
The primary output <sample_name>_convict_results.tsv provides the normalized coverage depth estimates in the last column. A companion script convict_concat.py can be used to concatenate multiple convict output files with default summary of the normalized coverage depths for each gene of interest. 

The variant identification component is still being tweaked as well as the ONT output parameterization using kmerresistance. 

The **tmp** directory has a number of alignment files as well as primary kmerresistance output to parse through for the user to troubleshoot or use to identify any potential gene variants further. 

If any comments, questions, or feedback please email me at wcshropshire@mdanderson.org
