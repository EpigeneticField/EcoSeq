# EcoSeq

Version 3.7

December 18, 2021

by Satoshi Yamashita (1) and Yuyu Liu (1)

Based on prior work by Naoko Iida (2)

(1) Division of Epigenomics, National Cancer Center Research Institute, Tokyo, Japan

(2) Division of Genome Analysis Platform Development, National Cancer Center Research Institute, Tokyo, Japan


1. Glossary

Eco-seq: Enzymatically Cleaved and Optimal Sequencing

UMI: Unique Molecular Identifier, random 12-mer in this study

OS: Original Sequence

Duplex Sequencing (https://github.com/Kennedy-Lab-UW/Duplex-Sequencing/blob/master/README.md)

 SSCS: Single Stranded Consensus Sequence

 DCS: Duplex Consensus Sequence

 Family: A group of reads sharing the same tag sequence


2. Dependencies

This bash script requires following programs with versions as specified:

Python >3.7 or >2.7

Duplex Sequencing software Version 3.0 July 11, 2016 Programs (our bash script has not yet adapted to the new Duplex Sequencing Pipeline software) https://github.com/Kennedy-Lab-UW/Duplex-Sequencing

R >3.6

gzip

fastp

bwa

samtools >1.9

popoolation2 https://github.com/popgenvienna/popoolation2

GenomonFisher https://github.com/Genomon-Project/GenomonFisher

Pysam (required for Duplex Sequencing software)

BioPython (required for Duplex Sequencing software)


3. Setup

Ecoseq_v37.sh

Make sure that the bash script is in the folder where you want your outputs.

Also make sure that the source folder (named "source") is in your working directory.

It should contain the following files:

EcoseqR_v21.R

dbSnp153Common3.bed

hg38.ecoseq.blacklist.v1.bed

Duplex Sequencing softwares

mpileup2sync.jar (from popoolation2)

hg38.fa


4. Input file preparation

Illumina 150bp paired-end reads:

Sample_a_1.fq.gz

Sample_a_2.fq.gz


Create an input file (txt format) containing the list of sample names (prefix identical to the fastq files). Leave a blank line at the end. 

Example of input file for the following fastq files: Sample_a_1.fq.gz, Sample_a_2.fq.gz, Sample_b_1.fq.gz, Sample_b_2.fq.gz, Sample_c_1.fq.gz, Sample_c_2.fq.gz:

Sample_a
Sample_b
Sample_c


5. Run the script as follows:

./Ecoseq_v37.sh <path to reference genome for bwa> <file with list of samples> <number of threads>

Example: ./Ecoseq_v37.sh /Volumes/Ueda01/2001_Epi_duplexSeq/bwa/hg38 samples.txt 8


6. Output file descriptions:

File Description                                |File name
------------------------------------------------|----------------------------------------------------------------------
Fastq files tagged:                             |Sample_a.seq1.smi.fq, Sample_a.seq2.smi.fq
Data files for UMI counting:                    |Sample_a.tagstats, Sample_a_data.txt, Sample_a.png, Sample_a.zoom.png
Fastp reports:                                  |Sample_a.f5.report.html, Sample_a.f5an0.uniq1.dcs.report.html
BAM file containing uniquely mapped reads:      |Sample_a.f5.uniq1.pe.sorted.bam
BAM files containing SSCS:                      |Sample_a.f5.uniq1.sscs_LCC.bam, Sample_a.f5.uniq1.sscs_NM.bam, Sample_a.f5.uniq1.sscs.bam, Sample_a.f5.uniq1.sscs.sorted.bam
Tagcounts file:                                 |Sample_a.f5.uniq1.pe.tagcounts
Tagstats file:                                  |Sample_a.f5.uniq1.pe.tagstats
BAM file containing DCS:                        |Sample_a.f5.uniq1.dcs.bam
Fastq files containing DCS:                     |Sample_a.f5.uniq1.dcs.r1.fq, Sample_a.f5.uniq1.dcs.r2.fq
Fastq files containing DCS trimmed:             |Sample_a.f5an0.uniq1.dcs.r1.fq, Sample_a.f5an0.uniq1.dcs.r2.fq
BAM file containing DCS sorted:                 |Sample_a.f5an0.uniq1.dcs.sorted.bam
BAM file containing DCS filtered:               |Sample_a.f5an0.uniq1.dcs.fil.bam
List of all variants:                           |Sample_a.f5an0.uniq1.dcs.fil.mutcount.txt
Bed file of all variants:                       |Sample_a.f5an0.uniq1.dcs.fil.mutcount.bed (use in R script)
Sequencing depth of DCS reads:                  |Sample_a.f5an0.uniq1.dcs.fil.mpileup.txt
Depth of DCS reads filtered by the blacklist:   |Sample_a.f5an0.uniqb.dcs.fil.mpileup.rdc.bed (use in R script)
Depth of OS reads filtered by the blacklist:    |Sample_a.f5an0.uniqb.dcs.fil.mpileup.pe.mpileup.txt (use in R script)
Bed files of analyzed regions with DCS depth(x):|Sample_a.uniqb.dcs.mpileup.x.20.bed
Mutation lists with DCS depth (x):              |Sample_a.uniqb.x.20.mu.txt
Mutation lists after removal of SNPs:           |Sample_a.uniqb.x.20.mu1.txt
Summary table for mutations:                    |Sample_a.uniqb.mu.results.txt
Summary table for sequencing depth:             |Sample_a.uniqb.depth.results.txt
