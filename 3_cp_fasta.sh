#!/bin/sh


source_dir="/xdisk/uschuch/corykeith/ViCAT_P003/out/oriented_alphas/partials"
NCBI_REF="/xdisk/uschuch/corykeith/BLAST/NCBI_virus/viraldb_1line.fsa"
while read line
    do
    DEST_DIR="/xdisk/uschuch/corykeith/ViCAT_P003/out/oriented_alphas/partials/correction/$line"
    echo $line
    echo $DEST_DIR
    for file in "$source_dir"/*.fasta; do
        if grep -q "$line" "$file"; then
            cp "$file" "$DEST_DIR"
            ACCESSION=${file#*NODE_*_}
            ACCESSION=${ACCESSION%_oriented.fasta}
            echo $ACCESSION >> "${DEST_DIR}/refseq_list.txt"
        fi
    done
    while IFS= read -r LINE
    do
    echo $LINE
    grep -A 1 $LINE $NCBI_REF > $DEST_DIR/${LINE}.fasta
    done < "${DEST_DIR}/refseq_list.txt"
done < sorted_alpha_samp_names.txt
