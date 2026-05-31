#!/bin/bash

# program to process fastq data with the goal
# (1) to select sequences with constant regions with some errors tolerated
# (2) to trimm the constant regions from the region of interest with some errors tolerated 
# (3) to derive the copy number of the region of interest and plot the result in a histogram

# Print the current working directory
echo "The current working directory is: $(pwd)"

Folder="data/FASTQ/ml_assessment_library"

for file in "$Folder"/*paired_ends_merged.fastq; do
    if [ -e "$file" ]; then
        base=$(basename "$file" _paired_ends_merged.fastq)
        echo "$base"
        # filter sequences with both left and right flanking regions with some errors tolerated and extract the variable region
        # changed min_overlap from 21 to 5 because I believe that a lot of longer variants are getting lost due to the processing
        cutadapt -g "TATCTCTCAAAGACTATTAAC;min_overlap=21;required...ACGCTAAAATTCAGTGTGGCC;min_overlap=21;required" -o data/FASTQ/${base}_filtered_trimmed.fastq --discard-untrimmed $file
    else
        echo "No .fastq files found in $Folder"
    fi
done
echo "finished processing"

    

