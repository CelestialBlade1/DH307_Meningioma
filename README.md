# DH307_Meningioma

# RNA Seq Pipeline Docs | Vijay
meningioma AND [RNA-Seq]

## Data Sample for Preliminary Set Up:
GSE136661: 
Obtaining Data Sample: http://www.ebi.ac.uk/ena/data/view/SRX6777597

## Anaconda Environment: pipeline
Source: https://docs.anaconda.com/free/anaconda/packages/using-r-language/

### Creating an environment with R
Download and install Anaconda.
Create a new conda environment with all the r-essentials conda packages built from CRAN:
conda create -n pipeline r-essentials r-base
Activate the environment:
conda activate pipeline
List the packages in the environment:
conda list
The list shows that the package r-base is installed and r is listed in the build string of the other R packages in the environment.
Anaconda Navigator, the Anaconda graphical package manager and application launcher, creates R environments by default.
FastQC
Install fastqc Using apt-get on Ubuntu
Update apt database with apt-get using the following command.
sudo apt-get update

After updating apt database, We can install fastqc using apt-get by running the following command:
sudo apt-get -y install fastqc

### Anaconda Installation of fastqc
Source: https://anaconda.org/bioconda/fastqc
conda install -c bioconda fastqc

### Trimmomatic
Conda Installation
Source: https://anaconda.org/bioconda/trimmomatic
conda install -c bioconda trimmomatic

### Hisat2
Conda Installation
Source: https://anaconda.org/bioconda/hisat2
conda install -c bioconda hisat2

### Stringtie
Conda Installation
Source: https://anaconda.org/bioconda/stringtie
conda install -c bioconda stringtie

### edgeR
Conda Installation
Source:https://anaconda.org/bioconda/bioconductor-edger
conda install -c bioconda bioconductor-edger

### Deseq2
Conda Installation (**)
Source: https://anaconda.org/bioconda/bioconductor-deseq2
conda install -c bioconda bioconductor-deseq2

### Samtools
Conda Installation
Source: https://anaconda.org/bioconda/samtools
conda install -c bioconda samtools

### Grch38 Reference Genome
Download and unzip reference: 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
Getting the annotated file: https://www.biostars.org/p/174331/ 
