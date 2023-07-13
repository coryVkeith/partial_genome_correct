#!/bin/sh

DIR="/xdisk/uschuch/corykeith/ViCAT_P003/out/oriented_alphas/partials/"
OUT_FILE="samp_names.txt"
for file in ${DIR}*.fasta
    do
    # Get the filename without the path
    filename=$(basename "$file")
    filename1="${filename%_NODE*}"
    # Write the filename to the output file
    echo "$filename1" >> "$OUT_FILE"
done

# De-duplicate the output file
sorted_file="sorted_alpha_samp_names.txt"
sort -u "$OUT_FILE" > "$sorted_file"

# Remove the original output file
#rm "$OUT_FILE"
