# Call pipeRunner.sh with correct arguments

import time
import sys
import argparse
import os
import numpy as np

parser = argparse.ArgumentParser(description='Run pipeline')
parser.add_argument('--input', type=str, help='Input text file containing list of SRA', required=True)
parser.add_argument('--output', type=str, help='Output directory', required=True)
parser.add_argument('--tempdir', type=str, help='Temporary directory', required=False, default='./tmp')
parser.add_argument('--genomedir', type=str, help='Directory containing genome fasta and annotation', required=False, default='./genome')
parser.add_argument('--adapterdir', type=str, help='Directory containing adapter fasta', required=False, default='./Trimmomatic-0.39/adapters')
parser.add_argument('--remove', type=str, help='Remove intermediate files', required=False, default='yes', choices=['yes', 'no'])
parser.add_argument('--progresslog', type=str, help='Progress log file', required=False, default='progress.log')

if __name__ == '__main__':
    args = parser.parse_args()
    print("Hello!!")
    REMOVE = True if args.remove == 'yes' else False
    FASTQ_DIR = "./fastq"
    # Create output directory
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    num_files = 0
    with open(args.input, 'r') as f:
        num_files = len([line for line in f])

    # Open input file and for each line, run the pipeline
    with open(args.input, 'r') as f:
        i = 0
        for line in f:
            # print("HERE??")
            i += 1
            # Create temporary directory
            if not os.path.exists(args.tempdir):
                # print("HERE??")
                os.makedirs(args.tempdir)
            # Create fastq directory
            if not os.path.exists(FASTQ_DIR):
                # print("HERO??")
                os.makedirs(FASTQ_DIR)
                    
            sra = line.strip()
            print("Processing SRA: " + sra)
            print(f"!!!!!!!! File {i}/{num_files} !!!!!!!!")
            # Run pipeline
            start = time.time()
            os.system("./pipeRunner.sh " + sra + " " + FASTQ_DIR + " " + args.output + " " + args.tempdir + " " + args.genomedir)
            end = time.time()
            print("Time taken: " + str(end - start) + " seconds")
            # Remove intermediate files
            if REMOVE:
                # Remove the entire temporary directory
                os.system("rm -rf " + args.tempdir)
                # Remove fastq files
                os.system("rm -rf " + FASTQ_DIR + "/*")
            
            # Write progress to log file
            with open(args.progresslog, 'a') as f:
                f.write(f"{i}/{num_files} {sra}\n")

