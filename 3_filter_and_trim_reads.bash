#!/bin/bash

# program to process fastq data with the goal
# (1) to select sequences with constant regions with some errors tolerated
# (2) to trimm the constant regions from the region of interest with some errors tolerated 

# Print the current working directory
echo "The current working directory is: $(pwd)"

Folder="$1"
OutputFolder="${2:-data/FASTQ}"

mkdir -p "$OutputFolder"

for file in "$Folder"/*paired_ends_merged.fastq; do
    if [ -e "$file" ]; then
        base=$(basename "$file" _paired_ends_merged.fastq)
        echo "$base"
        # filter sequences with both left and right flanking regions with some errors tolerated and extract the variable region
        cutadapt -g "TATCTCTCAAAGACTATTAAC;min_overlap=21;required...ACGCTAAAATTCAGTGTGGCC;min_overlap=21;required" -o "$OutputFolder/${base}_filtered_trimmed.fastq" --discard-untrimmed $file
    else
        echo "No .fastq files found in $Folder"
    fi
done
echo "finished processing"

    

