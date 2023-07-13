#!/bin/sh
set -u
#
# Checking args
#
source ./config.sh
#
#makes sure sample file is in the right place
#
if [[ ! -f "$PROFILE" ]]; then
    echo "$PROFILE does not exist. Please provide the path for a metagenome profile. Job terminated."
    exit 1
fi
#
# creates outdir
#
if [[ ! -d "$OUT_DIR" ]]; then
    echo "$OUT_DIR does not exist. The folder was created."
    mkdir $OUT_DIR
fi
## This builds individual directories from the list of samples in the profile list
while read SAMPLE; do
  echo "$SAMPLE"
    export SAMPLE_DIR="$OUT_DIR/$SAMPLE"
    export RAOUT="$OUT_DIR/relative_abundance"
    export BOWOUT="$SAMPLE_DIR/bowtie2/unused_reads"
    export ALNMNTOUT="$SAMPLE_DIR/bowtie2/alignments/stats"
    export CONSENSUS="$SAMPLE_DIR/bowtie2/consensus"
    export SPADESOUT="$SAMPLE_DIR/spades"
    export BLAST="$SAMPLE_DIR/blast"    
    export CONTIGFILT="$SAMPLE_DIR/spades_filtered"
    export ALLSPADES="$OUT_DIR/all_spades"
    export ALLBLAST="$OUT_DIR/all_blast"
    export CONTIGS="$OUT_DIR/virus_contigs"
    export BOWTIE2INDEX="$SAMPLE_DIR/blast_ref_fasta/bowtie2index"
    export ALLCONSENSUS="$OUT_DIR/all_consensus"
    export ERRORS="$OUT_DIR/errors"
    export COV="$SAMPLE_DIR/bowtie2/alignments/coverage/graphs"
    export ALLCOV="$OUT_DIR/all_cov"
#    init_dir "$SAMPLE_DIR" "$RAOUT" "$ALLCOV" "$BOWOUT"  "$ALNMNTOUT" "$SPADESOUT" "$BLAST" "$CONTIGFILT" "$ALLSPADES" "$ALLBLAST" "$CONTIGS" "$BOWTIE2INDEX" "$CONSENSUS" "$ALLCONSENSUS" "$ERRORS" "$COV" 
done <$PROFILE

# Job submission 
ARGS="-p $QUEUE --account=$GROUP --mail-user $MAIL_USER --mail-type $MAIL_TYPE"

# how many jobs in the array
export NUM_JOB=$(wc -l < "$PROFILE")

###
### SPAdes de novo and filters contigs
###
#Prog1="de_novo"
#export STDERR_DIR1="$SCRIPT_DIR/err/$Prog1"
#export STDOUT_DIR1="$SCRIPT_DIR/out/$Prog1"
#init_dir "$STDERR_DIR1" "$STDOUT_DIR1"
#
#echo "$SPADES"
#echo "launching $SCRIPT_DIR/run_spades.sh as a job."
#
######-a tells the number of jobs to submit to the PBS array
#JOB_ID=`sbatch $ARGS --export=ALL,SEQKIT=$SEQKIT,CSVTK=$CSVTK,AllSPADES=$ALLSPADES,OUT_DIR=$OUT_DIR,RAW=$RAW,STDERR_DIR1=$STDERR_DIR1,STDOUT_DIR1=$STDOUT_DIR1,PROFILE=$PROFILE,SPADES=$SPADES,SPADESOUT=$SPADESOUT,ERRORS=$ERRORS --job-name $Prog1 -o $STDOUT_DIR1/output.%a.out -e $STDERR_DIR1/err.%a.out -a 1-$NUM_JOB $SCRIPT_DIR/run_spades.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
####
#### blast filtered contigs
####
#Prog2="blast"
#export STDERR_DIR2="$SCRIPT_DIR/err/$Prog2"
#export STDOUT_DIR2="$SCRIPT_DIR/out/$Prog2"
#init_dir "$STDERR_DIR2" "$STDOUT_DIR2"
#
#echo " launching $SCRIPT_DIR/run_blastn.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,OUT_DIR=$OUT_DIR,ALLBLAST=$ALLBLAST,STDERR_DIR2=$STDERR_DIR2,STDOUT_DIR2=$STDOUT_DIR2,PROFILE=$PROFILE,BLASTDB=$BLASTDB,ERRORS=$ERRORS --job-name $Prog2 -o $STDOUT_DIR2/output.%a.out -e $STDERR_DIR2/err.%a.out --dependency=afterok:$PREV_JOB_ID -a 1-$NUM_JOB $SCRIPT_DIR/run_blastn.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
#####
##### This submits a code that writes accession numbers of blast hits to a list and does the same for contigs
#####
#Prog3="write_files"
#export STDERR_DIR3="$SCRIPT_DIR/err/$Prog3"
#export STDOUT_DIR3="$SCRIPT_DIR/out/$Prog3"
#init_dir "$STDERR_DIR3" "$STDOUT_DIR3"
#
#echo " launching $SCRIPT_DIR/run_collect_contigs_accessions.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,OUT_DIR=$OUT_DIR,STDERR_DIR3=$STDERR_DIR3,STDOUT_DIR3=$STDOUT_DIR3,PROFILE=$PROFILE,ERRORS=$ERRORS,BLAST1line=$BLAST1line --job-name $Prog3 -o $STDOUT_DIR3/output.%a.out -e $STDERR_DIR3/err.%a.out --dependency=afterok:$PREV_JOB_ID -a 1-$NUM_JOB $SCRIPT_DIR/run_collect_contigs_accessions.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
####
#### This submits a code that gathers the contigs from the filtered contigs list from those that have blast hits to a virus.
####
#Prog4="virus_contigs"
#export STDERR_DIR4="$SCRIPT_DIR/err/$Prog4"
#export STDOUT_DIR4="$SCRIPT_DIR/out/$Prog4"
#init_dir "$STDERR_DIR4" "$STDOUT_DIR4"
#
#echo " launching $SCRIPT_DIR/run_gather_contigs.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,SAMPLE_DIR=$SAMPLE_DIR,SCRIPT_DIR=$SCRIPT_DIR,OUT_DIR=$OUT_DIR,STDERR_DIR4=$STDERR_DIR4,STDOUT_DIR4=$STDOUT_DIR4,PROFILE=$PROFILE,GATHER=$GATHER,ERRORS=$ERRORS --job-name $Prog4 -o $STDOUT_DIR4/output.%a.out -e $STDERR_DIR4/err.%a.out --dependency=afterok:$PREV_JOB_ID -a 1-$NUM_JOB $SCRIPT_DIR/run_gather_contigs.sh`
##
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
####
#### Submits Code to remove redundancy in the fasta files used to build bowtie2 indexes for read mapping. The threshold is set to 88.5% as determined experimentally for begomoviruses.
####
#Prog5="CD-HIT"
#export STDERR_DIR5="$SCRIPT_DIR/err/$Prog5"
#export STDOUT_DIR5="$SCRIPT_DIR/out/$Prog5"
#init_dir "$STDERR_DIR5" "$STDOUT_DIR5"
#
#echo " launching $SCRIPT_DIR/run_cdhit.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,SAMPLE_DIR=$SAMPLE_DIR,OUT_DIR=$OUT_DIR,STDERR_DIR5=$STDERR_DIR5,STDOUT_DIR5=$STDOUT_DIR5,PROFILE=$PROFILE,GATHER=$GATHER,ERRORS=$ERRORS --job-name $Prog5 -o $STDOUT_DIR5/output.%a.out -e $STDERR_DIR5/err.%a.out --dependency=afterok:$PREV_JOB_ID -a 1-$NUM_JOB $SCRIPT_DIR/run_cdhit.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
####
#### This code builds a bowtie2 index from the combined fasta file of the accessions that the samples contigs hit and for the scaffolds themselves
####
#Prog6="bowtie2"
#export STDERR_DIR6="$SCRIPT_DIR/err/$Prog6"
#export STDOUT_DIR6="$SCRIPT_DIR/out/$Prog6"
#init_dir "$STDERR_DIR6" "$STDOUT_DIR6"
#
#echo " launching $SCRIPT_DIR/run_bowtie2buildnmap.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,BAMTOOLS=$BAMTOOLS,SAMPLE_DIR=$SAMPLE_DIR,ALLCONSENSUS=$ALLCONSENSUS,OUT_DIR=$OUT_DIR,STDERR_DIR6=$STDERR_DIR6,STDOUT_DIR6=$STDOUT_DIR6,PROFILE=$PROFILE,ALLBLAST=$ALLBLAST,BOWTIE=$BOWTIE,ERRORS=$ERRORS --job-name $Prog6 --dependency=afterok:$PREV_JOB_ID -o $STDOUT_DIR6/output.%a.out -e $STDERR_DIR6/err.%a.out -a 1-$NUM_JOB $SCRIPT_DIR/run_bowtie2buildnmap.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
####
#### This code splits the BAM file of the combined index to individual files and then calls stats
####
#Prog7="split"
#export STDERR_DIR7="$SCRIPT_DIR/err/$Prog7"
#export STDOUT_DIR7="$SCRIPT_DIR/out/$Prog7"
#init_dir "$STDERR_DIR7" "$STDOUT_DIR7"
#
#echo " launching $SCRIPT_DIR/run_split.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,SAMPLE_DIR=$SAMPLE_DIR,BAMTOOLS=$BAMTOOLS,OUT_DIR=$OUT_DIR,STDERR_DIR7=$STDERR_DIR7,STDOUT_DIR7=$STDOUT_DIR7,PROFILE=$PROFILE,ALLBLAST=$ALLBLAST,BOWTIE=$BOWTIE,ERRORS=$ERRORS --job-name $Prog7 --dependency=afterok:$PREV_JOB_ID -o $STDOUT_DIR7/output.%a.out -e $STDERR_DIR7/err.%a.out -a 1-$NUM_JOB $SCRIPT_DIR/run_split.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
####
#### This code writes the coverage information to graphs for each taxa in the sample.
####
#Prog8="cov"
#export STDERR_DIR8="$SCRIPT_DIR/err/$Prog8"
#export STDOUT_DIR8="$SCRIPT_DIR/out/$Prog8"
#init_dir "$STDERR_DIR8" "$STDOUT_DIR8"
#
#echo " launching $SCRIPT_DIR/run_split.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,SAMPLE_DIR=$SAMPLE_DIR,SCRIPT_DIR=$SCRIPT_DIR,BAMTOOLS=$BAMTOOLS,OUT_DIR=$OUT_DIR,STDERR_DIR8=$STDERR_DIR8,STDOUT_DIR8=$STDOUT_DIR8,PROFILE=$PROFILE,ALLBLAST=$ALLBLAST,BOWTIE=$BOWTIE,ERRORS=$ERRORS --job-name $Prog8 --dependency=afterok:$PREV_JOB_ID -o $STDOUT_DIR8/output.%a.out -e $STDERR_DIR8/err.%a.out -a 1-$NUM_JOB $SCRIPT_DIR/run_covgraph_bed.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
####
#### This code writes the coverage information to graphs for each taxa in the sample.
####
#Prog9="taxmat"
#export STDERR_DIR9="$SCRIPT_DIR/err/$Prog9"
#export STDOUT_DIR9="$SCRIPT_DIR/out/$Prog9"
#init_dir "$STDERR_DIR9" "$STDOUT_DIR9"
#
#echo " launching $SCRIPT_DIR/run_split.sh in queue"
#echo "previous job ID $PREV_JOB_ID"
#
#JOB_ID=`sbatch $ARGS --export=ALL,SAMPLE_DIR=$SAMPLE_DIR,SCRIPT_DIR=$SCRIPT_DIR,BAMTOOLS=$BAMTOOLS,OUT_DIR=$OUT_DIR,STDERR_DIR9=$STDERR_DIR9,STDOUT_DIR9=$STDOUT_DIR9,PROFILE=$PROFILE,ALLBLAST=$ALLBLAST,BOWTIE=$BOWTIE,ERRORS=$ERRORS --job-name $Prog9 --dependency=afterok:$PREV_JOB_ID -o $STDOUT_DIR9/output.%a.out -e $STDERR_DIR9/err.%a.out -a 1-$NUM_JOB $SCRIPT_DIR/run_taxmat.sh`
#
#if [ "${JOB_ID}x" != "x" ]; then
#        JOB_ID=${JOB_ID#"Submitted batch job "}
#
#        echo Job: \"$JOB_ID\"
#        PREV_JOB_ID=$JOB_ID
#else
#        echo Problem submitting job. Job terminated.
#        exit 1
#fi
#
###
### This code writes the coverage information to graphs for each taxa in the sample.
###
Prog10="partial_correct"
export STDERR_DIR10="$SCRIPT_DIR/err/$Prog10"
export STDOUT_DIR10="$SCRIPT_DIR/out/$Prog10"
init_dir "$STDERR_DIR10" "$STDOUT_DIR10"

echo " launching $SCRIPT_DIR/run_split.sh in queue"
#echo "previous job ID $PREV_JOB_ID"

JOB_ID=`sbatch $ARGS --export=ALL,MINIMAP=$MINIMAP,SAMPLE_DIR=$SAMPLE_DIR,SCRIPT_DIR=$SCRIPT_DIR,OUT_DIR=$OUT_DIR,STDERR_DIR10=$STDERR_DIR10,STDOUT_DIR10=$STDOUT_DIR10,PROFILE=$PROFILE,ALLBLAST=$ALLBLAST,BOWTIE=$BOWTIE,ERRORS=$ERRORS --job-name $Prog10 -o $STDOUT_DIR10/output.%a.out -e $STDERR_DIR10/err.%a.out -a 1-$NUM_JOB $SCRIPT_DIR/run_alpha_correct.sh`

if [ "${JOB_ID}x" != "x" ]; then
        JOB_ID=${JOB_ID#"Submitted batch job "}

        echo Job: \"$JOB_ID\"
        PREV_JOB_ID=$JOB_ID
else
        echo Problem submitting job. Job terminated.
        exit 1
fi

