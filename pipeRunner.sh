#!/bin/bash


# Pipeline: FastQC -> Trimmomatic -> hisat2 -> StringTie -> DESeq2/edgeR
# The pipeline gets executed for only one sample at a time
# See pipeline.py for the code that calls this script


SRR=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
TEMPFILES_DIR=$4
GENOME_DIR=$5
ADAPTER_DIR=$6

# Time the pipeline
start_time=$(date +%s)


echo "#############################################################"
echo "Pipeline executing on $SRR"

echo "#############################################################"
echo "Running FastQC: Step 1/4"
fasterq-dump --split-files $SRR -O $INPUT_DIR -p -v

echo "#############################################################"
echo "Running Trimmomatic: Step 2/4"
trimmomatic PE -threads 12 $INPUT_DIR/$SRR"_1.fastq" $INPUT_DIR/$SRR"_2.fastq" -baseout $TEMPFILES_DIR/$SRR".fastq.gz" ILLUMINACLIP:$ADAPTER_DIR/TruSeq3-PE.fa:4:30:10 MINLEN:36 TRAILING:3 LEADING:3

echo "#############################################################"
echo "Running hisat2: Step 3/4"
hisat2 -x $GENOME_DIR/genome -1 $TEMPFILES_DIR/$SRR"_1P.fastq.gz" -2 $TEMPFILES_DIR/$SRR"_2P.fastq.gz" -S $TEMPFILES_DIR/$SRR.sam -p 12

echo "#############################################################"
echo "Running StringTie: Step 4/4"
samtools sort -o $TEMPFILES_DIR/$SRR.bam $TEMPFILES_DIR/$SRR.sam
samtools index $TEMPFILES_DIR/$SRR.bam
mkdir -p $OUTPUT_DIR/$SRR
stringtie $TEMPFILES_DIR/$SRR.bam -o $OUTPUT_DIR/$SRR/$SRR.gtf -p 8 -l $SRR -G $GENOME_DIR/genome_ann.gtf -e

end_time=$(date +%s)

# Print the time taken to run the pipeline
echo "Time taken to run pipeline upto stringtie: $((end_time-start_time)) seconds"
echo "Pipeline complete for $SRR"
echo "#############################################################"