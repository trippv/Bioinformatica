#!/bin/bash
#############################################

# Counts RNAseq fragments (read pairs) mao to each annotated gene in the genome
# We  need to provide our GTF formatted annotation so reads can be assigned to gene features.

# Directivas
#SBATCH --job-name=06_HTcount
#SBATCH --output=06_HTcount-%j.log # agregar num de trabajo
#SBATCH --error=06_HTcount-%j.err
#SBATCH --nodes 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00 #numero maximo
#SBATCH -p cicese #cola

#################################################################
# Generate Counts 
#################################################################

# htseq runs in conda 3.5 
module load python-3.5-anaconda 

#source /LUSTRE/apps/Anaconda/2019/miniconda3/etc/profile.d/conda.sh

#### paths to prorams
SAMTOOLS=/LUSTRE/apps/bioinformatica/samtools-1.7/bin
HTSEQ=/home/rgomez/bin

## path to reads
INDIR=../../../Trimmed_READS/Trimmed_nonRNA_P

# path to genome
GENOME=../genome/Crassostrea_gigas.GCA902806645v1.dna.toplevel.fa

#path to annotation
GFF=../genome/Crassostrea_gigas.GCA902806645v1.54.gff3

#Path to bam alignments
ALIGNMENTS=../alignments

#Counts directory
OUTDIR=../counts
mkdir -p $OUTDIR


### Alternative run:
#htseq-count -f bam --stranded=no GCF_902806645.1_cgigas_uk_roslin_v1_genomic.bam GCF_902806645.1_cgigas_uk_roslin_v1_genomic.gff > couts.txt


cd ${SLURM_SUBMIT_DIR}

#module load perl-5.30.0

echo "Using bam alignmend from $ALIGMENTS"


#Get the sample name from original reads
for inFile1 in $INDIR/*_R1_p.fq; do
   inFile2=${inFile1/_R1_p.fq/_R2_p.fq}
   outBase=${inFile1%%_R1_p.fq}

    #Extract sample name
   file=${inFile2##*/}
   SAMPLE=${file%%_R*_p.fq}

  echo "Using file ${SAMPLE}.bam"


#Runhtseq count

$HTSEQ/htseq-count -f bam \
    --stranded=no \
    -r pos \
    -t exon \
    --idattr=exon  \
    -f bam $ALIGNMENTS/${SAMPLE}.bam \
    $GFF \
    > $OUTDIR/${SAMPLE}.counts


done
#Now we will use the program htseq-count to count how many RNA fragments (i.e. read pairs) map to each annotated gene in the 
#genome. Again, we'll use our accession list and the program parallel. We also need to provide our GTF formatted annotation so 
#reads can be assigned to gene features.

#-s no indicates we're using an unstranded RNA-seq library.
#-r pos tells htseq-count that our BAM file is coordinate sorted.
#-f bam indicates that our input file is in BAM format.
#-t feature type (3rd column in ghe gff file) to be sued, all features of other type are ignored (default exon)


## Run success
exit
