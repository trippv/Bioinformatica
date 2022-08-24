#!/bin/bash
#############################################
# Script to get the reference genome and annotation

# Directivas
#SBATCH --job-name=get-genome
#SBATCH --output=get-genome-%j.log # agregar num de trabajo
#SBATCH --error=get-genome-%j.err
#SBATCH --nodes 1
#SBATCH --mem=20GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00 #numero maximo
#SBATCH -p cicese #cola


# Print history to log file

echo 'download reference genome'

date

echo "download genome:"
echo "http://ftp.ensemblgenomes.org/pub/metazoa/release-54/fasta/crassostrea_gigas/dna/Crassostrea_gigas.GCA902806645v1.dna.toplevel.fa.gz"

#################################################################
# Download genome and annotation from ENSEMBL
#################################################################

# load software
# module load samtools/1.12

# output directory
GENOMEDIR=../genome
mkdir -p $GENOMEDIR

# we're using the crassotrea gigas genome GCA902806645v1 (UK Roslin)


    # we'll download the genome, GTF annotation and transcript fasta
    # https://useast.ensembl.org/Fundulus_heteroclitus/Info/Index

    # Source: Carolina Peñaloza[1], The Roslin Institute and Royal (Dick) School of Veterinary Studies, The University of Edinburgh, 2022

# download the genome
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-54/fasta/crassostrea_gigas/dna/Crassostrea_gigas.GCA902806645v1.dna.toplevel.fa.gz
# decompress it
gunzip Crassostrea_gigas.GCA902806645v1.dna.toplevel.fa.gz

# download the GTF annotation
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-54/gff3/crassostrea_gigas/Crassostrea_gigas.GCA902806645v1.54.gff3.gz
# decompress it
gunzip Crassostrea_gigas.GCA902806645v1.54.gff3.gz



# generate simple samtools fai indexes 
samtools faidx Crassostrea_gigas.GCA902806645v1.dna.toplevel.fa.gz


# move everything to the genome directory
mv Crassotrea* $GENOMEDIR