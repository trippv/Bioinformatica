#!/bin/bash
#############################################
# Script to make the genome index with HISAT

# Directivas
#SBATCH --job-name=03_HISATndex
#SBATCH --output=03_HISATindex-%j.log # agregar num de trabajo
#SBATCH --error=03_HISATindex-%j.err
#SBATCH --nodes 1
#SBATCH --mem=20GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00 #numero maximo
#SBATCH -p cicese #cola


date

#################################################################
# Index the Genome
#################################################################


HTSEQ=/LUSTRE/apps/bioinformatica/hisat2-2.1.0


## Usage: hisat2-build [options]* <reference_in> <ht2_base>

## <reference_in> Reference sequences to be aligned to
## <ht2_base> basename of the index files 

cd ${SLURM_SUBMIT_DIR}


# load python environment

module load python-3.7-anaconda 


# input/output directories
OUTDIR=../genome/hisat2_index
mkdir -p $OUTDIR

GENOME=../genome/Crassostrea_gigas.GCA902806645v1.dna.toplevel.fa

hisat2-build -p 16 $GENOME $OUTDIR/Cgig