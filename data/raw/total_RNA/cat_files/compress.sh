#!/bin/bash

# Author: Jorge A. Castro
# Version: 2.0.0
# File compression after sample concatenation.

directory=/Users/Parsimony/GitHub_Projects/mosquito_virome_yucatan_LEVE/data/raw/total_RNA/cat_files

for file in "$directory"/*fastq; do
	if [ -f "$file" ]; then
		echo "Compressing: $file" 	
		gzip "$file"
	fi
done

echo "Process completed..."
