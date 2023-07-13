#!/bin/sh

# Your job will use 1 node, 28 cores, and 168gb of memory total.
#SBATCH -N 1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=6G
#SBATCH -t 18:00:00



## calls the profile variable to pull sample names from a list iteratively
echo $SLURM_ARRAY_TASK_ID >> taskid_check.txt
export SMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $PROFILE | tail -n 1`
echo $SMPLE >> smple_check.txt
#file="${RAW}/*.gz"
for file in ${RAW}/*.gz; do
    if echo $file | grep -q "$SMPLE"; then
       filename=$(basename "$file")
       BASE="${filename%_R*_001.fastq.gz}"
       F1="${BASE}_R1_001.fastq.gz"
       R1="${BASE}_R2_001.fastq.gz"
       echo $F1 $R1 >> read_check.txt
    fi
done

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate ViCAT

HELPER_DIR="$OUT_DIR/oriented_helpers/partials/correction"
BETA_DIR="$OUT_DIR/oriented_betas/partials/correction"
ALPHA_DIR="$OUT_DIR/oriented_alphas/partials/correction"
#mkdir -p $HELPER_DIR
#mkdir -p $BETA_DIR
mkdir -p $ALPHA_DIR
mkdir -p "$OUT_DIR/oriented_helpers/partials/correction/${SMPLE}/bowtie2/idx"
mkdir -p "$OUT_DIR/oriented_helpers/partials/correction/combined_fasta"
mkdir -p "$OUT_DIR/oriented_helpers/partials/correction/combined_cov"
#mkdir -p "$OUT_DIR/oriented_betas/partials/correction/${SMPLE}/bowtie2/idx"
#mkdir -p "$OUT_DIR/oriented_betas/partials/correction/combined_fasta"
#mkdir -p "$OUT_DIR/oriented_betas/partials/correction/combined_cov"
#mkdir -p "$OUT_DIR/oriented_alphas/partials/correction/${SMPLE}/bowtie2/idx"
#mkdir -p "$OUT_DIR/oriented_alphas/partials/correction/combined_fasta"
#mkdir -p "$OUT_DIR/oriented_alphas/partials/correction/combined_cov"

## Counting number of references in folder
REF_COUNT=0
for f in ${ALPHA_DIR}/${SMPLE}/NC_*.fasta
     do REF_COUNT=$(( REF_COUNT + 1 ))
     echo $f >> ${ALPHA_DIR}/${SMPLE}/REF_NAMES.txt
done
CONTIG_NUM=$(( REF_COUNT + 1 ))

## Building References Array for mauve with multiple references.
REFS=()

while IFS= read -r current_line; do
    REFS+=("$current_line")
done < ${ALPHA_DIR}/${SMPLE}/REF_NAMES.txt

###Combine partial fasta files into one multiple fasta
cat ${ALPHA_DIR}/${SMPLE}/*NODE* > ${ALPHA_DIR}/${SMPLE}/${SMPLE}_combined_partials.fasta


$MINIMAP -k 10 -w 5 -B 0 -ado ${REFS[0]} ${ALPHA_DIR}/${SMPLE}/${SMPLE}_combined_partials.fasta > ${ALPHA_DIR}/${SMPLE}/${SMPLE}_minimap.sam
     # Convert sam to bam
     samtools view -S -b ${ALPHA_DIR}/${SMPLE}/${SMPLE}_minimap.sam > ${ALPHA_DIR}/${SMPLE}/${SMPLE}_minimap.bam
     # Sort the alignment
     samtools sort ${ALPHA_DIR}/${SMPLE}/${SMPLE}_minimap.bam -o ${ALPHA_DIR}/${SMPLE}/${SMPLE}_minimap_sorted.bam
     # Get consensus fastq file
     samtools mpileup -B -uf ${REFS[0]} ${ALPHA_DIR}/${SMPLE}/${SMPLE}_minimap_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${ALPHA_DIR}/${SMPLE}/${SMPLE}_tmp_consensus.fastq
     #change header to name of contig instead of reference
     seq_ref1=$(basename "${REFS[0]}")
     seq_ref1=${seq_ref1%.fasta}
     SEQ_DESC1="${SMPLE}_$seq_ref1"
     echo $SEQ_DESC1 $seq_ref1 >> ${ALPHA_DIR}/${SMPLE}/seq_desc.txt
     cat ${ALPHA_DIR}/${SMPLE}/${SMPLE}_tmp_consensus.fastq | seqkit replace -p $seq_ref1 -r $SEQ_DESC1 > ${ALPHA_DIR}/${SMPLE}/${SMPLE}_consensus.fastq
 #    rm ${BAM1}_tmp_consensus.fastq
     # Convert .fastq to .fasta
     seqtk seq -a ${ALPHA_DIR}/${SMPLE}/${SMPLE}_consensus.fastq > ${ALPHA_DIR}/${SMPLE}/${SMPLE}_${seq_ref1}_partial_corrected.fasta


###bowtie2 build index of filtered references above
if [[ -f "$ALPHA_DIR/${SMPLE}/${SMPLE}_${seq_ref1}_partial_corrected.fasta" ]]; then
    Combined_partials="$ALPHA_DIR/${SMPLE}/${SMPLE}_${seq_ref1}_partial_corrected.fasta"
    $BOWTIE/bowtie2-build -f $Combined_partials $ALPHA_DIR/${SMPLE}/bowtie2/idx/${SMPLE}_idx
    else
    echo "Error in run_bowtie2buildnmap.sh on partial genome reconstruction, bowtie2. No combined partials fasta file. Stopping at build." >> $ERRORS/${SMPLE}_error.log
fi
##bowtie2 reference mapping to indexes
if [[ -f "$ALPHA_DIR/${SMPLE}/bowtie2/idx/${SMPLE}_idx.1.bt2" ]]; then    
    $BOWTIE/bowtie2  -x $ALPHA_DIR/${SMPLE}/bowtie2/idx/"${SMPLE}_idx" -q -1 ${RAW}/$F1 -q -2 ${RAW}/$R1 -S $ALPHA_DIR/${SMPLE}/bowtie2/${SMPLE}.SAM
    else
    echo "Error in run_bowtie2buildnmap.sh, bowtie2. No index. Stopping at mapping" >> $ERRORS/${SMPLE}_error.log
fi
#Samtools convert to sorted bam, extract consensus, and remove sam
if [[ -f "${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.SAM" ]]; then
    samtools view -S -b ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.SAM > ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.bam
    samtools sort ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.bam  -o ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}_sorted.bam
    samtools mpileup --max-depth 0 -uf $Combined_partials ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.fastq
    seqtk seq -a ${BETA_DIR}/${SMPLE}/bowtie2/${SMPLE}.fastq > ${BETA_DIR}/${SMPLE}/bowtie2/${SMPLE}_partial_corrected.fasta

#_____________________________________________________________________________________________________________________________________|    
###    This code will copy the consensus sequences to a common folder and remove the .sam and unsorted .bam files to help save space. |
###                                     Comment this code if you want to keep intermediate files.                                     |
###                                                                                                                                   |
#_____________________________________________________________________________________________________________________________________|    
    cp ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}_partial_corrected.fasta "${ALPHA_DIR}/combined_fasta"
    rm ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.SAM
    rm ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.bam
    else
    echo "Error in run_bowtie2buildnmap.sh, samtools. No sam alignment file found. Stopping at bam file conversion." >> $ERRORS/${SMPLE}_error.log
fi

if [[ -f "${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}_sorted.bam" ]]; then
        f="${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}_sorted.bam"
        bedtools genomecov -bg -ibam $f > ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.coverage
        sed -i -r 's/(\s+)?\S+//2' ${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.coverage
    touch ${ALPHA_DIR}/${SMPLE}/bowtie2/coverage.txt
else
    echo "Error in calculating coverage. No reference sorted bam files found." >> $ERRORS/${SMPLE}_error.log
fi

PLOT_NAME="${SMPLE}_partial_corrected"

if [[ -f "${ALPHA_DIR}/${SMPLE}/bowtie2/${SMPLE}.coverage" ]]; then
    module load R
    Rscript $SCRIPT_DIR/cov_for_graph.R $ALPHA_DIR $SMPLE --save
    touch ${ALPHA_DIR}/${SMPLE}/bowtie2/svg.txt
else
    echo "Error in generating coverage graphs. No coverage information used to build graph" >> $ERRORS/${SMPLE}_error.log
fi

if [[ -f "${ALPHA_DIR}/${SMPLE}/bowtie2/svg.txt" ]]; then
    cp ${ALPHA_DIR}/${SMPLE}/bowtie2/*.svg ${ALPHA_DIR}/combined_cov
else
    echo "Error in copying coverage graphs to common folder. No coverage graph." >> $ERRORS/${SMPLE}_error.log
fi
