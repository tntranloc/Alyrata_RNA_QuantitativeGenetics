#!/bin/bash

# Define input and output file
input_file="variants_table.tsv"
output_file="pi_ratio_output.tsv"

# Print header to output file
echo -e "Gene\tGeneId\tTranscriptId\tπi\tπs\tπi/πs" > $output_file

# Skip the first line (header) and process each row
tail -n +2 $input_file | while IFS=$'\t' read -r gene gene_id transcript_id biotype high low moderate modifier \
        utr3 premature_start utr5 downstream initiator intron missense noncoding splice_acceptor \
        splice_donor splice_region start_lost stop_gained stop_lost stop_retained synonymous upstream; do

    # Extract πi (missense variants) and πs (synonymous variants)
    pi_i=$missense
    pi_s=$synonymous

    # Calculate πi/πs, handle division by zero
    if [ "$pi_s" -ne 0 ]; then
        pi_ratio=$(echo "scale=4; $pi_i / $pi_s" | bc)
    else
        pi_ratio="undefined"
    fi

    # Append results to the output file
    echo -e "$gene\t$gene_id\t$transcript_id\t$pi_i\t$pi_s\t$pi_ratio" >> $output_file
done

# Inform the user
echo "Calculation complete! Results saved in $output_file"
