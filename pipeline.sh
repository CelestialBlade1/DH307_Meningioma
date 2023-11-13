# Pipeline: FastQC -> Trimmomatic -> hisat2 -> StringTie -> DESeq2/edgeR
# This script does not take command line arguments. All parameters are set in the script.
# Please read the documentation for directory structure and naming conventions.
# To run: bash pipeline.sh and give permission to execute if necessary (chmod +x pipeline.sh)
# We only emulate the first 3 steps of the pipeline for now.


INPUT_DIR="./data"
FASTQC_DIR="./fastqc"
TRIM_DIR="./trim"

# Directory for hg38 reference genome in .fa format
GENOME_DIR="./genome"

echo "This bash acript assumes paired end data with _1.fastq.gz and _2.fastq.gz"
echo "as input in the INPUT_DIR. Please ensure that or you will face errors"

sleep 3
 
# Find number of files in INPUT_DIR
num_files=$(ls -1 $INPUT_DIR/*_1.fastq.gz | wc -l)

# Initialize counter
count=0

# The input directory contains Paired End (PE) reads with names <NAME>_1.fastq.gz and <NAME>_2.fastq.gz
# We want to extract tha <NAME> part of the file name to use as the base name for all output files

declare -A unique_names
count=0
# Iterate over the files in the directory
for file in "$INPUT_DIR"/*_1.fastq.gz; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        
        # Extract the base name of the file (without extension)
        base_name=$(basename "$file" | cut -d '_' -f 1)

        # Add the base name to an associative array to keep track of unique names
        unique_names["$base_name"]="$base_name"

        # Debug: Print the file being processed and the extracted base name
        echo "Processing file: $file, Base Name: $base_name"
        count=$((count+1))
    fi
done

# Debug: Print the list of unique names
# echo "Unique Names:"
# for count in "${!unique_names[@]}"; do
#     echo "${unique_names[$count]}"
# done

# Time the pipeline
start_time=$(date +%s)

# echo "Files to process: $files"
count=0
for file in "${!unique_names[@]}"; do
    # Increment counter
    count=$((count+1))
    # file=${unique_names[$count]}

    # Print start of pipeline step
    echo "#############################################################"
    echo "Pipeline step $count of $num_files: Executing on $file"
    echo "#############################################################"

    # FastQC step
    fastqc $INPUT_DIR/$file"_1.fastq.gz" $INPUT_DIR/$file"_2.fastq.gz" -o $FASTQC_DIR
    echo "----------FastQC step complete----------"
    
    # Run Trimmomatic
    echo "----------Running Trimmomatic----------"
    mkdir $TRIM_DIR/$file
    trimmomatic PE -threads 4 $INPUT_DIR/$file"_1.fastq.gz" $INPUT_DIR/$file"_2.fastq.gz" -baseout $TRIM_DIR/$file/$file.fastq.gz ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE.fa:4:30:10 MINLEN:36 TRAILING:3 LEADING:3
    echo "----------Trimmomatic step complete----------"

    # Run hisat2
    echo "----------Running hisat2----------"
    mkdir -p ./hisat2/$file
    # Hisat on paired end reads
    hisat2 -x $GENOME_DIR/genome -1 $TRIM_DIR/$file/$file"_1P.fastq.gz" -2 $TRIM_DIR/$file/$file"_2P.fastq.gz" -S ./hisat2/$file/$file.sam
    # Hisat on single end reads
    # hisat2 -x $GENOME_DIR/genome -U $TRIM_DIR/$file/$file_1U.fastq.gz,$TRIM_DIR/$file/$file_2U.fastq.gz -S ./hisat2/$file/$file.sam  # Currently not working
    echo "----------hisat2 step complete----------"

    # Run samtools
    echo "----------Running samtools----------"
    samtools sort -o ./hisat2/$file/$file.bam ./hisat2/$file/$file.sam
    samtools index ./hisat2/$file/$file.bam
    echo "----------samtools step complete----------"

    # Run StringTie
    echo "----------Running StringTie----------"
    mkdir -p ./stringtie/$file
    stringtie ./hisat2/$file/$file.bam -o ./stringtie/$file/$file.gtf -p 8 -l $file
    echo "----------StringTie step complete----------"


done

# End timing the pipeline
end_time=$(date +%s)

# Print the time taken to run the pipeline
echo "Time taken to run pipeline upto stringtie: $((end_time-start_time)) seconds"

# Run prepDE.py
echo "----------Running prepDE.py----------"
mkdor ./prepDE
python3 prepDE.py -i ./stringtie -g ./prepDE/gene_count_matrix.csv -t ./prepDE/transcript_count_matrix.csv
echo "----------prepDE.py step complete----------"

# Run DESeq2.R
echo "----------Running DESeq2.R----------"
Rscript DESeq2.R
echo "----------DESeq2.R step complete----------"

# Another endtime
end_time2=$(date +%s)

# Print the time taken to run the pipeline
echo "Time taken to run pipeline upto DESeq2: $((end_time2-start_time)) seconds"

