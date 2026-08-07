#!/bin/bash

Folder="$1"
OutputFolder="$2"

if [ -z "$OutputFolder" ]; then
    echo "Usage: $0 <input_folder> <output_folder>"
    exit 1
fi

# Populate the array with file groups
for dir in "$Folder"/*; do
    if [ -d "$dir" ]; then
        declare -A fileGroups  # Create an associative array to hold file groups for each directory
        dirName=$(basename "$dir")
        outDir="$OutputFolder/$dirName"
        mkdir -p "$outDir"
        for file in "$dir"/*.fq; do
            base=$(basename "$file")
            # Remove L01, L02, and L03 from the base name
            base=${base/_L01_/_}
            base=${base/_L02_/_}
            base=${base/_L03_/_}
            group=${base%%.*}
            fileGroups["$group"]+="$file "
        done
        # Concatenate the files in each group
        for group in "${!fileGroups[@]}"; do
            echo "Creating ${outDir}/${dirName}_${group}.fq"
            cat ${fileGroups["$group"]} > "${outDir}/${dirName}_${group}.fq"
        done
        unset fileGroups  # Unset the fileGroups array
    else
        echo "No directories found in $Folder"
    fi
done
