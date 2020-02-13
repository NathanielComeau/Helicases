# Helicases

```txt
`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
   `=`,'=/     `=`,'=/     `=`,'=/     `=`,'=/
     y==/        y==/        y==/        y==/
   ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.
,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
```

# BLOCKAGES: Where we need help from Ibrahim

- *Nat: Sample data. I grabbed some human genome data from https://www.internationalgenome.org/data-portal/sample/HG00114 ; is this appropriate?*
- Nat: Our paper is to be in the Nature Bioinformatics format right?
- Is our project looking for dead cameras, etc and compressing based on that? Re-ordering quality scores, or a combination of the two? **Big Picture Help.**
- More of the same question: Is this problem just re-ordering the quality scores in the best possible way? Is this a separate problem than looking for broken cameras, etc and compressing based on that? Nat: can Ibrahim compare and contrast, in Scalce, the compression method used for reads and the one used for quality scores? *Is our project essentially trying to adapt the loss-less method used for reads to quality scores, using flow-cell etc. info?*
- How does Scalce define an optimal ordering?
- Nat: I'm finding a hard time figuring out where exactly compression is happening in the Scalce source code (for reads).
- Nat: FastQC? Anything useful in there?

### Answers:

- *Nat: Show Ibrahim full_genome_exploration/side_by_side.txt . What are some potential clustering methods we can use on this?*
    - **Answer: k-means Clustering.**

## Immediate To-Dos

- Better tool to extract identifiers and quality scores from a giant FASTQ file.

- K-Means Clustering

- Do a quick sanity check on identifier components to explore their properties.
  1) Is name & instrument the same for all identifiers?
  2) Are all their lengths equal?

## Eventual To-Dos

### C++ Implementation

We need a working idea first, but eventually this will need implemented at production scale. The C++ source code for scalce is quite nice but
it has lots of extra features we might not need. I'm wondering what can we take from Scalce for our prototype?

### Not important- Nat tried to install Ibrahim's Benchmarking Toolkit 

Tried to follow the instructions https://github.com/sfu-compbio/compression-benchmark/blob/master/aws.md

, got as far as vagrant up before erroring out because of an "image not found" error.

## Random Things Nat Learned

- It turns out read lengths are not always 100 bp long. Who would have thunk it!

