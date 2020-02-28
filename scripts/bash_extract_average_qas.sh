#!/bin/bash

# Extracts average quality scores from each line of a file of quality scores. 

# Check that a filename was given
if [ -n "$1" ];
then
    : #proceed as normal
else
    echo "Usage: $0 <quality score file>"
    exit 1
fi


awk 'BEGIN {
    FS=""
    for(n=0;n<256;n++)ord[sprintf("%c",n)]=n
}
{
	count = 0;
    for (i=1;i<=NF;i++) count += ord[$i];
    count /= NF;
    count -= 33;
    print count;

}

END {
}' $1