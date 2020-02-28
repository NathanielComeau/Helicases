#!/bin/bash

# Extracts median quality scores from each line of a file of quality scores. 

# Check that a filename was given
if [ -n "$1" ];
then
    : #proceed as normal
else
    echo "Usage: $0 <quality score file>"
    exit 1
fi


gawk 'BEGIN {
    FS=""
    for(n=0;n<256;n++)ord[sprintf("%c",n)]=n;

    QA_OFFSET = 33;
}
{
    for (i=1;i<=NF;i++) ar[i] = ord[$i];
    n = asort(ar);
    print (ar[int(n/2)] - QA_OFFSET);
}

END {
}' $1