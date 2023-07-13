#!/bin/sh


while read line;
    do
    mkdir -p /xdisk/uschuch/corykeith/ViCAT_P003/out/oriented_alphas/partials/correction/$line
done < sorted_alpha_samp_names.txt


