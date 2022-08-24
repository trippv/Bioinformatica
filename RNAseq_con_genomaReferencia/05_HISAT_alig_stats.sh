#!/bin/bash
#############################################
#Map RNAseq read to genome with HiSat

#Once we have created the index, then we align the reads to the reference genome with HISAT2 

# Directivas
#SBATCH --job-name=05_alignQC
#SBATCH --output=05_alignQC-%j.log # agregar num de trabajo
#SBATCH --error=05_alignQC-%j.err
#SBATCH --nodes 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00 #numero maximo
#SBATCH -p cicese #cola


##################################
# calculate stats on alignments
##################################

# calculate alignment statistics for each bam file using samtools

# load software--------------------------------------------------------------------------
### path to programs ########
SAMTOOLS=/LUSTRE/apps/bioinformatica/samtools-1.7/bin
HISAT=/LUSTRE/apps/bioinformatica/hisat2-2.1.0
#############################################

# input, output directories--------------------------------------------------------------

INDIR=../alignments
OUTDIR=../samtools_stats
MULTIQC=../samtools_stats_multiqc
mkdir -p $MULTIQC
mkdir -p $OUTDIR


# samtools bam statistics----------------------------------------------------------------

$SAMTOOLS/samtools stats $INDIR/Ch14_1.bam >$OUTDIR/Ch14_1.stats
$SAMTOOLS/samtools stats $INDIR/Ch14_2.bam >$OUTDIR/Ch14_2.stats
$SAMTOOLS/samtools stats $INDIR/Ch14_3.bam >$OUTDIR/Ch14_3.stats
$SAMTOOLS/samtools stats $INDIR/EP08_1.bam >$OUTDIR/EP08_1.stats
$SAMTOOLS/samtools stats $INDIR/EP08_2.bam >$OUTDIR/EP08_2.stats
$SAMTOOLS/samtools stats $INDIR/EP08_3.bam >$OUTDIR/EP08_3.stats
$SAMTOOLS/samtools stats $INDIR/Se17_1.bam >$OUTDIR/Se17_1.stats
$SAMTOOLS/samtools stats $INDIR/Se17_2.bam >$OUTDIR/Se17_2.stats
$SAMTOOLS/samtools stats $INDIR/Se17_3.bam >$OUTDIR/Se17_3.stats

# put the basic stats all in one file.---------------------------------------------------

# bash array containing all stats file names
FILES=($(find $OUTDIR -name "*.stats" | sort))

# initialize big table for all samples with row names
grep "^SN" ${FILES[0]} | cut -f 2 > $OUTDIR/SN.txt

# loop over each stats file, add data column to main table
for file in ${FILES[@]}
do paste $OUTDIR/SN.txt <(grep ^SN $file | cut -f 3) > $OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt
	echo $file
done

# add a header with sample names
cat \
<(echo ${FILES[@]} | sed 's,samtools_stats/,,g' | sed 's/.stats//g' | sed 's/ /\t/g') \
$OUTDIR/SN.txt \
>$OUTDIR/SN2.txt && \
	mv $OUTDIR/SN2.txt $OUTDIR/SN.txt



#### make multiqc repor

module load python-2.7-anaconda 
source activate multiqc_py2.7

multiqc -f -o $MULTIQC $OUTDIR