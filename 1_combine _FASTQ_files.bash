#!/bin/bash

Folder="$1"

# Populate the array with file groups
for dir in "$Folder"/*; do
    if [ -d "$dir" ]; then
        declare -A fileGroups  # Create an associative array to hold file groups for each directory
        dirName=$(basename "$dir")
        for file in "$dir"/*.fq; do
            base=$(basename "$file")
            # Remove L01, L02, and L03 from the base name
            base=${base/*_L01_*-/}
            base=${base/*_L02_*-/}
            base=${base/*_L03_*-/}
            group=${base%%.*}
            fileGroups["$group"]+="$file "
        done
        # Concatenate the files in each group
        for group in "${!fileGroups[@]}"; do
            echo "Creating ${dir}/${dirName}_${group}.fq"
            cat ${fileGroups["$group"]} > "${dir}/${dirName}_${group}.fq"
        done
        unset fileGroups  # Unset the fileGroups array
    else 
        echo "No directories found in $Folder"
    fi
done
