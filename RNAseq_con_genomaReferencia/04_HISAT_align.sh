#!/bin/bash
#############################################
#Map RNAseq read to genome with HiSat

#Once we have created the index, then we align the reads to the reference genome with HISAT2 

# Directivas
#SBATCH --job-name=04_HISATalign
#SBATCH --output=04_HISATalign-%j.log # agregar num de trabajo
#SBATCH --error=04_HISATalign-%j.err
#SBATCH --nodes 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00 #numero maximo
#SBATCH -p cicese #cola

#################################################################
# Align reads to genome
#################################################################

INDIR=../../../Trimmed_READS/Trimmed_nonRNA_P #trimmed reads path. 
OUTDIR=../alignments #Output path
genome=../genome
mkdir -p $OUTDIR

INDEX=../genome/hisat2_index/Cgig #path to index
EXONS=../genome/exonsFile.table

### path to programs ########
SAMTOOLS=/LUSTRE/apps/bioinformatica/samtools-1.7/bin
HISAT=/LUSTRE/apps/bioinformatica/hisat2-2.1.0
#############################################

## Usage: hisat2-build [options]* <reference_in> <ht2_base>

## <reference_in> Reference sequences to be aligned to
## <ht2_base> basename of the index files 


############### make exons splice site ########################
# The python command sends an error so this can be done manually#
echo "make exon file"
cat $genome/Crassostrea_gigas.GCA902806645v1.54.gff3 | awk '{if ($3=="exon") {print $1"\t"$4-1"\t"$5-1}}'  > $genome/exonsFile.table
#################################################################


### Alternative run:
#srun hisat2 -x index/Cgig_Roslin_index --known-splicesite-infile exonsFile.table -p 24 -1 ../../All_paired_R1.fq.gz -2 ../../All_paired_R2.fq.gz 2> Hisat_mapping.err 1> Hisat_map

cd ${SLURM_SUBMIT_DIR}


# load python environment
module load py
module load perl-5.30.0

# print info to log file
echo "Using reads from $INDIR"
echo "in file $OUTDIR"

for inFile1 in $INDIR/*_R1_p.fq; do
   inFile2=${inFile1/_R1_p.fq/_R2_p.fq}

   echo "mapping $inFile1 + $inFile2"
   outBase=${inFile1%%_R1_p.fq}
   #outBase=${outBase/$INDIR/$OUTDIR}

    #Extract sample name
   file=${inFile2##*/}
   SAMPLE=${file%%_R*_p.fq}

   echo $SAMPLE


#Run Hisat mapping to reads
$HISAT/hisat2 -x $INDEX \
--known-splicesite-infile $EXONS \
-p 24 \
-1 $INDIR/${SAMPLE}_R1_p.fq   \
-2 $INDIR/${SAMPLE}_R2_p.fq | \
 
 $SAMTOOLS/samtools view -@ 8 -S -h -u - | \
 $SAMTOOLS/samtools sort -@ 8 -T $SAMPLE - >$OUTDIR/$SAMPLE.bam


# index bam files

$SAMTOOLS/samtools index $OUTDIR/$SAMPLE.bam 

done
#The | is the pipe. It tells linux to use the output of the command to the left of the pipe as the input for the command 
#to the right. You can chain together many commands using pipes. samtools view converts the SAM file produced by hisat2 to 
#uncompressed BAM. -S indicates the input is SAM format. -h indicates the SAM header should be written to the output. 
#-u indicates that uncompressed BAM format should be written (no need to compress until the end). 
#- indicates samtools should take input from the pipe. samtools sort sorts and compressed the file. 
#Finally, the -@ option can be used to allocate additional threads to be used for compression,
#-T gives a temporary file prefix.


## Run success
exit


