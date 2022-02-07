##EcoSeq_official_v37##
##./EcoSeq_official_v37.sh /Users/epigenome/hg38/hg38.fa sample2.txt 4

#!/bin/bash

I_DIR=$(pwd)
S_DIR=${I_DIR}/source
REF=$1
TH=$3

if [ ! -d ${S_DIR} ] ; then
	echo -e "ERROR: Cannot find source files! Please check your working space...\n       Exiting program..."
	exit 1
fi

while read line; do
	NAME=${I_DIR}/${line}
	echo -e "$(date +"%T") Analyzing $line...\n         Decompressing fastq.gz files...\n"
	# Decompress
	if [ ! -f ${NAME}_1.fq.gz ] || [ ! -f ${NAME}_2.fq.gz ] ; then
		echo -e "WARNING: No such file: ${line}! \n         Please check your working space...\n         Skipping this sample..."
		break
	fi
	gzip -dc ${NAME}_1.fq.gz > ${NAME}_R1.fastq
	gzip -dc ${NAME}_2.fq.gz > ${NAME}_R2.fastq

	# Identify barcode
	echo -e "$(date +"%T") Identifying barcode...\n"
	if [ ! -f ${S_DIR}/tag_to_header.py ] ; then
		echo  -e "ERROR: Cannot find tag_to_header.py script! Please check your source folder...\n  Exiting program..."
		exit 1
	fi
	python ${S_DIR}/tag_to_header.py --infile1 ${NAME}_R1.fastq --infile2 ${NAME}_R2.fastq --outprefix ${line} --tagstats --spacerlen 5 --taglen 12 --filtspacer 'TGACT'

	# Fastp QC & trim
	echo -e "\n$(date +"%T") QC & trimming...\n"
	if ! type "fastp" > /dev/null ; then
		echo -e "ERROR: fastp is not installed! Please install all dependencies required...\n       Exiting program..."
		exit 1
	fi
	fastp -i ${NAME}.seq1.smi.fq -I ${NAME}.seq2.smi.fq -o ${NAME}.f5.seq1.smi.fq -O ${NAME}.f5.seq2.smi.fq -h ${NAME}.f5.report.html -f 5 -F 5 -A -w $TH

	# BWA alignment
	echo -e '\n'$(date +"%T")' Mapping to genome...\nThis might take a while...\n'
	if ! type "bwa" > /dev/null ; then
		echo -e "ERROR: bwa is not installed! Please install all dependencies required...\n       Exiting program..."
		exit 1
	fi
	bwa mem -t $TH $REF ${NAME}.f5.seq1.smi.fq ${NAME}.f5.seq2.smi.fq > ${NAME}.f5.pe.sam
	echo -e '\n'$(date +"%T")' Cleaning mapped results...\n'
	grep -v -E -e '\bXA:Z:' -e '\bSA:Z:' ${NAME}.f5.pe.sam > ${NAME}.f5.uni.pe.sam
	if ! type "samtools" > /dev/null ; then
		echo -e "ERROR: samtools is not installed! Please install all dependencies required...\n       Exiting program..."
		exit 1
	fi
	samtools view -@ $TH -b -q 1 ${NAME}.f5.uni.pe.sam > ${NAME}.f5.uniq1.pe.bam
	samtools sort -@ $TH -m 4G -O bam -o ${NAME}.f5.uniq1.pe.sorted.bam ${NAME}.f5.uniq1.pe.bam
	samtools index -@ $TH ${NAME}.f5.uniq1.pe.sorted.bam

	# Make SSCS
	echo -e '\n'$(date +"%T")' Making SSCS...\n'
	if [ ! -f ${S_DIR}/ConsensusMaker.py ] ; then
		echo -e "ERROR: Cannot find ConsensusMaker.py script! Please check your source folder...\n       Exiting program..."
		exit 1
	fi
	python ${S_DIR}/ConsensusMaker.py --infile ${NAME}.f5.uniq1.pe.sorted.bam --tag_file ${NAME}.f5.uniq1.pe.tagcounts --tag_stats ${NAME}.f5.uniq1.pe.tagstats --outfile ${NAME}.f5.uniq1.sscs.bam --minmem 3 --read_length 128 --cut_off 0.7 --Ncut_off 0.3 --read_type dp --filt n
	samtools sort -@ $TH -m 4G -O bam -o ${NAME}.f5.uniq1.sscs.sorted.bam ${NAME}.f5.uniq1.sscs.bam
	if [ ! -f ${NAME}.f5.uniq1.sscs.sorted.bam ] ; then
		echo -e "ERROR: Failed making consensus sequences for ${line}! \n         Skipping this sample..."
		break
	fi
	samtools index -@ $TH ${NAME}.f5.uniq1.sscs.sorted.bam

	# Make DCS
	echo -e '\n'$(date +"%T")' Making DCS...\n'
	if [ ! -f ${S_DIR}/DuplexMaker.py ] ; then
		echo -e 'ERROR: Cannot find DuplexMaker.py script! Please check your source folder...\n       Exiting program...'
		exit 1
	fi
	python ${S_DIR}/DuplexMaker.py --infile ${NAME}.f5.uniq1.sscs.sorted.bam --outfile ${NAME}.f5.uniq1.dcs.bam --Ncutoff 0.3 --readlength 128
	fastp -i ${NAME}.f5.uniq1.dcs.r1.fq -I ${NAME}.f5.uniq1.dcs.r2.fq -o ${NAME}.f5an0.uniq1.dcs.r1.fq -O ${NAME}.f5an0.uniq1.dcs.r2.fq -h ${NAME}.f5an0.uniq1.dcs.report.html -a GGATCAGTCA -w $TH -n 0
	bwa mem -t $TH $REF ${NAME}.f5an0.uniq1.dcs.r1.fq ${NAME}.f5an0.uniq1.dcs.r2.fq > ${NAME}.f5an0.uniq1.dcs.aln.sam
	samtools sort -@ $TH -m 4G -O bam -o ${NAME}.f5an0.uniq1.dcs.sorted.bam ${NAME}.f5an0.uniq1.dcs.aln.sam
	samtools index -@ $TH ${NAME}.f5an0.uniq1.dcs.sorted.bam
	
	# Filter DCS (properly paired)
	samtools view -@ $TH -f 2 -F 2048 -b ${NAME}.f5an0.uniq1.dcs.sorted.bam > ${NAME}.f5an0.uniq1.dcs.fil.bam
	samtools index -@ $TH ${NAME}.f5an0.uniq1.dcs.fil.bam

	# Detect variants
	echo -e $(date +"%T")' Detecting variants...\n'
	if ! type "fisher" > /dev/null ; then
		echo -e "ERROR: GenomonFisher is not installed! Please install all dependencies required...\n       Exiting program..."
		exit 1
	fi
	fisher single -1 ${NAME}.f5an0.uniq1.dcs.fil.bam -o ${NAME}.f5an0.uniq1.dcs.fil.mutcount.txt -s samtools -r ${S_DIR}/hg38.fa -m 0.0  -p 0.0 -d 1 -v 1 -e -S "-BQ 0"
	sed -e '1d' ${NAME}.f5an0.uniq1.dcs.fil.mutcount.txt | awk '{print $1,$2-1,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > ${NAME}.f5an0.uniq1.dcs.fil.mutcount.bed

	# Analyze and filter variants
	echo -e '\n'$(date +"%T")' Analyzing and filtering variants...\n'
	samtools mpileup ${NAME}.f5an0.uniq1.dcs.fil.bam > ${NAME}.f5an0.uniq1.dcs.fil.mpileup
	if ! type "java" > /dev/null ; then
		echo -e "ERROR: JAVA is not installed! Please install all dependencies required...\n       Exiting program..."
		exit 1
	fi
	if [ ! -f ${S_DIR}/mpileup2sync.jar ] ; then
		echo -e "ERROR: Cannot find mpileup2sync.jar script! Please check your source folder...\n       Exiting program..."
		exit 1
	fi
	java -ea -Xmx7g -jar ${S_DIR}/mpileup2sync.jar --input ${NAME}.f5an0.uniq1.dcs.fil.mpileup --output ${NAME}.f5an0.uniq1.dcs.fil.mpileup.txt --fastq-type sanger --threads $TH
	awk -F "[\t:]" 'BEGIN {OFS="\t"}{print $1,$2-1,$2,$4,$5,$6,$7,$8,$9}' ${NAME}.f5an0.uniq1.dcs.fil.mpileup.txt > ${NAME}.f5an0.uniq1.dcs.fil.mpileup.bed
	if ! type "bedtools" > /dev/null ; then
		echo -e "ERROR: bedtools is not installed! Please install all dependencies required...\n       Exiting program..."
		exit 1
	fi
	if [ ! -f ${S_DIR}/hg38.ecoseq.blacklist.v1.bed ] ; then
		echo -e "ERROR: Cannot find hg38.ecoseq.blacklist.v1.bed file! Please check your source folder...\n       Exiting program..."
		exit 1
	fi
	bedtools intersect -a ${NAME}.f5an0.uniq1.dcs.fil.mpileup.bed -b ${S_DIR}/hg38.ecoseq.blacklist.v1.bed -v > ${NAME}.f5an0.uniqb.dcs.fil.mpileup.rdc.bed
	samtools mpileup -a -l ${NAME}.f5an0.uniqb.dcs.fil.mpileup.rdc.bed ${NAME}.f5.uniq1.pe.sorted.bam > ${NAME}.f5an0.uniqb.dcs.fil.mpileup.pe.mpileup
	java -ea -Xmx7g -jar ${S_DIR}/mpileup2sync.jar --input ${NAME}.f5an0.uniqb.dcs.fil.mpileup.pe.mpileup --output ${NAME}.f5an0.uniqb.dcs.fil.mpileup.pe.mpileup.txt --fastq-type sanger --threads $TH
	rm ${NAME}_R1.fastq ${NAME}_R2.fastq ${NAME}.f5.seq1.smi.fq ${NAME}.f5.seq2.smi.fq ${NAME}.f5.pe.sam ${NAME}.f5.uni.pe.sam ${NAME}.f5.uniq1.pe.bam ${NAME}.f5an0.uniq1.dcs.aln.sam ${NAME}.f5an0.uniq1.dcs.fil.mpileup ${NAME}.f5an0.uniq1.dcs.fil.mpileup.bed ${NAME}.f5an0.uniqb.dcs.fil.mpileup.pe.mpileup
	echo -e '\n'$(date +"%T")' DONE for '$line'...\n'
	
	# Calculate through R script
	echo -e $(date +"%T")' Calculating using R script...\n'
	if ! type "R" > /dev/null ; then
		echo -e "ERROR: R is not installed! Please install all dependencies required...\n       Exiting program..."
		exit 1
	fi
	Rscript ${S_DIR}/EcoseqR_v21.R ${I_DIR} $line
	
	# Movie files to sample folder
	if [ ! -d ${I_DIR}/${line} ] ; then
		mkdir ${I_DIR}/${line}
	fi
	O_DIR=${I_DIR}/${line}
	mv ${I_DIR}/${line}*.* ${O_DIR}
done < "$2"