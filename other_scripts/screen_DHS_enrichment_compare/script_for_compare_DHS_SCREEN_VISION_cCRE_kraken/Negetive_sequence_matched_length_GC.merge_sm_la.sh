#!/bin/bash
sm_file=$1
la_file=$2
output_file=$3

paste $sm_file $la_file | awk -F '\t' -v OFS='\t' '{if (($8-$4)^2 <= ($16-$12)^2) print $1,$2,$3,$4,$5,$6,$7,$8; else print $9,$10,$11,$12,$13,$14,$15,$16}' > $output_file


