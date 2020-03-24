#!/bin/bash

# Parses a full genome.

# REQUIRES IN CURRENT DIRECTORY:
# - parseidentifiers.sh
# - bash_extract_average_qas.sh


# Check that a filename was given
if [ -n "$1" ];
then
    : #proceed as normal
else
    echo "Usage: $0 <fastq file>"
    exit 1
fi


# Optionally unzip the file
# gunzip $1

# Extract identifiers and quality scores
awk 'NR % 4 == 0' $1 > quality_scores.txt
awk 'NR % 4 == 1' $1 > identifiers.txt

# Parse identifiers into constituent components
./parseidentifiers.sh identifiers.txt
# Parse quality scores into average quality scores
./bash_extract_average_qas.sh quality_scores.txt > average_qas.txt

# MANUAL FROM HERE

# Bin data based on distribution of x and y
# Plot each bin and stitch them together into one whole human genome plot







