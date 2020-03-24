# Helicases

```txt
`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
   `=`,'=/     `=`,'=/     `=`,'=/     `=`,'=/    `=`,'=/     `=`,'=/     `=`,'=/     `=`,'=/
     y==/        y==/        y==/        y==/       y==/        y==/        y==/        y==/
   ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.   ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,-<=`.
,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
```

# BLOCKAGES: Where we need help from Ibrahim

- Nat: I found a conference proceeding link claiming to do something similar to what we're trying to do, but I can't actually find the conference proceedings anywhere. Can someone help me find them?
- Nat, maybe for Ashley? **Do we know which identifier field is which in the sample data Ibrahim gave us? Can Ashley double check the format of identifier fields?** Might be possible to figure this out by reading the "recommended-fastq" paper from 482c. We could maybe ask Ibrahim as well.
- Nat: can we learn something from previous work in quality score recalibration? Maybe there are some good insights there that we can use in making our model. This would involve looking up papers on quality score recalibration and asking Ibrahim if he has some guidance for where to start.
- Nat: How does Samcomp's model work?
	- *Figured it out*: See samcomp directory
- Nat: Our paper is to be in the Nature Bioinformatics format right?
	- Yup, see announcement on connex
- Is our project looking for dead cameras, etc and compressing based on that? Re-ordering quality scores, or a combination of the two? **Big Picture Help.**
- More of the same question: Is this problem just re-ordering the quality scores in the best possible way? Is this a separate problem than looking for broken cameras, etc and compressing based on that? Nat: can Ibrahim compare and contrast, in Scalce, the compression method used for reads and the one used for quality scores? *Is our project essentially trying to adapt the loss-less method used for reads to quality scores, using flow-cell etc. info?*
- How does Scalce define an optimal ordering?
- Nat: I'm finding a hard time figuring out where exactly compression is happening in the Scalce source code (for reads).
- Nat: FastQC? Anything useful in there?
- Nat: Should we use Ibrahim's benchmarking toolkit for this?

### Answers:

- *Nat: Show Ibrahim full_genome_exploration/side_by_side.txt . What are some potential clustering methods we can use on this?*
    - **Answer: k-means Clustering.**

## Immediate To-Dos

- Implementation of basic compression with a simple model.
- Implementation of compression with samcomp's model (see files in samcomp directory)

## Eventual To-Dos

- Different model for each x coordinate
- Different model for each x and y bin

- Nat: K-Means Clustering

- From Ibrahim: "arith. coding or rANS, or samcomp/O2 models"

### C++ Implementation

## Random Things Nat Learned

- It turns out read lengths are not always 100 bp long. Who would have thunk it!

