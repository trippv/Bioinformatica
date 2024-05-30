# Script para generar indice y alinear lecturas de RNAseq a genoma de referencia

## Descripción
* genoma de Haliotis fulgens


```
#!/bin/bash
#############################################
# Script to make the genome index with HISAT

# Directivas
#SBATCH --job-name=Hisat
#SBATCH --output=Hisat_align_Hruf-%j.log # agregar num de trabajo

#SBATCH --nodes 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00 #numero maximo
#SBATCH -p cicese #cola



## script for HIsat2 alignments to the Haliotis fulgens genome


########################################################################################################

# path to software
HISAT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/hisat2-2.2.1
STRINGTIE=/LUSTRE/bioinformatica_data/genomica_funcional/apps/stringtie
GFFCOMPARE=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffcompare
GFFREAD=/LUSTRE/bioinformatica_data/genomica_funcional/apps/gffread #para convertif gff a gtf
SAMTOOLS=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3


############## parameters which can be modifed ################

# set the master root folder for the analysis
ROOT=/LUSTRE/bioinformatica_data/genomica_funcional/tripp/2024_ROBERTO/COMPARATIVO_ALIGNMENTS

#path to fasq files
INDIR=$ROOT/ALLREADS

# Path to agenome
GENOME_PATH=/LUSTRE/bioinformatica_data/genomica_funcional/tripp/GENOMES/Hful
GENOME=$GENOME_PATH/SoftmaskedFilteredHalful_medaka.FINAL.fasta
ANNOTATION_GFF=$GENOME_PATH/FulgensAnnotation.gff3
ANNOTATION_GTF=$GENOME_PATH/FulgensAnnotation.gtf
HISAT_STATS=$ROOT/HISAT_STATS



# set input FastQ patterns
FQ_pattern_R1='R1_p.fq.gz'
FQ_pattern_R2='R2_p.fq.gz'

#################################################################################
# index folder
index=$ROOT/INDEX
index_name="Hful"

#################################################################################

### 1 Preparar el indice ----------------------------------------------
# En caso de que no se tenga un archivo de anotacon GTF y solo se tenga un GFF, se puede convertir usando gffread de cufflink.
# Esta es una situación comun en especial con genomas de NCBI

### 1.1 Convertir el archivo gff a gtf
echo "convert gff to gtf with gffred"

$GFFREAD/gffread $ANNOTATION_GTF -T -o $ANNOTATION_GTF


## 1.2 generar archivo de exones y splice
$HISAT/extract_splice_sites.py $ANNOTATION_GTF > $GENOME_PATH/splice_sites.txt
$HISAT/extract_exons.py $ANNOTATION_GTF > $GENOME_PATH/exons.txt


# This method provides the sames genes and coords as the pyton script but with no direction (+, -)
#echo "make exon file manually"
#cat $ANNOT | awk '{if ($3=="exon") {print $1"\t"$4-1"\t"$5-1}}'  > $exonFile

## 1.3 generar el indice con la informacion de splice y exonces
#echo "make index"
echo "generar indice"



# Build Bowtie2 index if not already present
if [ ! -f "$index/$index_name.1.bt2" ]; then
    $HISAT/hisat2-build -p 24 $GENOME $index/$index_name \
    --exon $GENOME_PATH/exons.txt \
    --ss $GENOME_PATH/splice_sites.txt 2> HISAT_align_Hful_index_build.err

else
    echo "index '$index_name' already exists."
fi




echo "parte II alinear lecturas a genoma"

### PArt 2. align reads to genome (index) ----------------------------------------------------------------
# print info to log file
echo "Using reads from $INDIR"
echo "in file $HISAT_STATS"


for inFile1 in $INDIR/*$FQ_pattern_R1; do
   inFile2=${inFile1/_$FQ_pattern_R1/_$FQ_pattern_R2}

   #inFile1 -> nombre completo del arhicvo F
   #inFile2 -> nombre completo del archivo R

   echo "mapping $inFile1 + $inFile2"
   outBase=${inFile1%%$FQ_pattern_R1}
   #outBase=${outBase/$INDIR/$OUTDIR}

    #Extract sample name
   file=${inFile1##*/}
   SAMPLE=${file%%_$FQ_pattern_R1}


   #######################################


#Run Hisat mapping to reads
$HISAT/hisat2 -x $index/$index_name \
-p 24 \
--dta \
--max-intronlen 20000 \
-1 $inFile1   \
-2 $inFile2 2> $HISAT_STATS/${SAMPLE}_align_stats.txt | \
$SAMTOOLS/samtools view -@ 8 -S -h -u - | \
$SAMTOOLS/samtools sort -@ 8 -T $HISAT_STATS/${SAMPLE} - > $HISAT_STATS/$SAMPLE.bam


# index bam files

$SAMTOOLS/samtools index $HISAT_STATS/$SAMPLE.bam



### Make stats of the aligment
$SAMTOOLS/samtools stats $HISAT_STATS/$SAMPLE.bam > $HISAT_STATS/$SAMPLE.stats

done






```

