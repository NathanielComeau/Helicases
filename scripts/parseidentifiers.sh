#!/bin/bash

# Check that a filename was given
if [ -n "$1" ];
then
    echo "Parsing $1"
else
    echo "Usage: $0 <identifier file"
    exit 1
fi

# If directory identifierComponents exists, remove it and its contents
if [ -d "identifierComponents" ]; then
  rm -r "identifierComponents"
fi
mkdir "identifierComponents"

# Extract each identifier field to a separate text file
echo "Starting read names"
awk -F'[ .:#/]' '{print $1}' $1 > identifierComponents/read_names.txt
echo "Starting read numbers"
awk -F'[ .:#/]' '{print $2}' $1 > identifierComponents/read_numbers.txt
echo "Starting instrument names"
awk -F'[ .:#/]' '{print $3}' $1 > identifierComponents/instrument_names.txt
echo "Starting flow cell lanes"
awk -F'[ .:#/]' '{print $4}' $1 > identifierComponents/flow_cell_lanes.txt
echo "Starting tile numbers"
awk -F'[ .:#/]' '{print $5}' $1 > identifierComponents/tile_numbers.txt
echo "Starting x_coords"
awk -F'[ .:#/]' '{print $6}' $1 > identifierComponents/x_coords.txt
echo "Starting y_coords"
awk -F'[ .:#/]' '{print $7}' $1 > identifierComponents/y_coords.txt
echo "Done!"
