# Helicases

# Nat's work on weekend of Feb 8th

Summary: I started to explore quality score data from a sample fastq file. Try running "cd quality_scores; python3 qa_manipulation.py" to see some pretty plots of quality scores from 1, 10, and 100 reads.


- Downloaded scalce and built it from source (didn't build on MacOS, had to build on Linux)
- Downloaded an e-coli genome from https://www.ncbi.nlm.nih.gov/assembly/GCF_000008865.2
- Whoops, figured out that was an assembly file not a file with a ton of reads in it...
- Downloaded the SRR2040271_1 and SRR2040271_2 files instead. Ran scalce on the first file and got a compression ratio of 7.16.

```txt
nat@nat-MacBookPro:~/Desktop/Helicases/testFASTQfiles$ scalce SRR2040271_1.fastq -o compressed.fastq
SCALCE 2.8 [pthreads; available cores=4]
** pigz not found in PATH. Using zlib compression instead. **
Buffer size: 131072K, bucket storage size: 4194304K
Allocating pool... OK!
Preprocessing FASTQ files ...
    Paired end #1, quality offset: 33
                   read length: 100
OK
Using 3 threads...
    Done with file SRR2040271_1.fastq, 13516445 reads found
    Created temp file 0
Merging results ... 1
Generating 1 scalce file(s) from temp #0 compression: gzip
shrink factor for 0: 1
Written 1524012 cores
    Read bit size:    paired end 0 = 0.99
    Quality bit size: paired end 0 = 2.40
Cleaning ...
Statistics:
    Total number of reads: 13516445
    Read length: first end 100
    Unbucketed reads count: 619, bucketed percentage 100.00
    Lossy percentage: 0
    Time elapsed:     00:02:46
    Compression time: 00:00:55
    Original size: 4318.45M, new size: 603.13M, compression factor: 7.16
Done!

```

- Tried to compress both paired end reads at once, got a compression ratio of 7.9:

```txt
nat@nat-MacBookPro:~/Desktop/Helicases/testFASTQfiles$ scalce SRR2040271_1.fastq -r -o paired_end_compressed -n my_library
SCALCE 2.8 [pthreads; available cores=4]
** pigz not found in PATH. Using zlib compression instead. **
Buffer size: 131072K, bucket storage size: 4194304K
Allocating pool... OK!
Preprocessing FASTQ files ...
    Paired end #1, quality offset: 33
                   read length: 100
    Paired end #2, quality offset: 33
                   read length: 100
OK
Using 3 threads...
    Done with file SRR2040271_1.fastq, 13516445 reads found
              file SRR2040271_2.fastq, 13516445 reads found
    Created temp file 0
Merging results ... 1
Generating 1 scalce file(s) from temp #0 using paired files, compression: gzip
shrink factor for 0: 1
Written 1524013 cores
    Read bit size:    paired end 0 = 0.99
    Quality bit size: paired end 0 = 2.40
shrink factor for 1: 1
Written 1524013 cores
    Read bit size:    paired end 1 = 1.26
    Quality bit size: paired end 1 = 2.10
Cleaning ...
Statistics:
    Total number of reads: 13516445
    Read length: first end 100
                 second end 100
    Unbucketed reads count: 619, bucketed percentage 100.00
    Lossy percentage: 0
    Time elapsed:     00:03:11
    Compression time: 00:01:25
    Original size: 8636.91M, new size: 1086.89M, compression factor: 7.95
Done!
```

- Read the scalce(http://bioinformatics.oxfordjournals.org/content/28/23/3051) paper for a while. Learned:
	1) Ibrahim probably built it!
	2) They implemented a basic lossy compression system for quality scores.


## I started exploring the quality score data

- Started with SRR2040271_1.fastq
- Extracted 1000 records with seqkit.
- Extracted quality scores (every 4th line) with awk (quality_scores.txt).
- Also made a file with all quality scores (all_quality_scores.txt)
- Plotted them with ./quap (used plotSingleQualityScore("quality_scores.txt", 0, 33))
- Saved results in a jupyter notebook (exploring_quality_scores.ipynb)


```bash
seqkit head -n 1000 SRR2040271_1.fastq > thousand_records.fastq
awk 'NR % 4 == 0' thousand_records.fastq > quality_scores.txt
awk 'NR % 4 == 0' SRR2040271_1.fastq > all_quality_scores.txt

```

## I also started exploring the identifier data (flow cell, etc.)

- Started with SRR2040271_1.fastq
- Extracted 1000 records with seqkit.
- Extracted names (every 4th line after first) with awk (quality_scores.txt).

```bash
seqkit head -n 1000 SRR2040271_1.fastq > thousand_records.fastq
awk 'NR % 4 == 1' thousand_records.fastq > identifiers.txt

``` 

**Shoot, it didn't work for this genome! The identifiers were totally different
then the [illumina sequence information](https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers):**

```bash
w134-87-145-178:flow_cell_lanes nat$ head identifiers.txt 
@SRR2040271.1 SN603_WBP007_8_1101_63.30_99.90 length=100
@SRR2040271.2 SN603_WBP007_8_1101_79.90_99.30 length=100
@SRR2040271.3 SN603_WBP007_8_1101_218.30_99.00 length=100
@SRR2040271.4 SN603_WBP007_8_1101_247.90_99.30 length=100
@SRR2040271.5 SN603_WBP007_8_1101_275.90_99.10 length=100
@SRR2040271.6 SN603_WBP007_8_1101_399.20_99.50 length=100
@SRR2040271.7 SN603_WBP007_8_1101_419.60_99.20 length=100
@SRR2040271.8 SN603_WBP007_8_1101_421.50_99.60 length=100
@SRR2040271.9 SN603_WBP007_8_1101_459.10_99.70 length=100
@SRR2040271.10 SN603_WBP007_8_1101_487.80_99.60 length=100
```

Where from the FASTQ wikipedia page we expect something like:

    - @HWUSI-EAS100R:6:73:941:1973#0/1
    - @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

From the wiki page identifiers take one of the following formats:

```text
<Instrument Name>:<tile number>:<x-coordinate>:<y-coordinate>...#<index number>/<member of pair>
```
Where everything after ... is complicateded, or:

```text
<Instrument Name>:<Run ID>:<Flowcell ID>:<Flowcell Lane>:<Tile Number>:<X Coord>:<Y Coord> ...
```
Again where everything after ... is complicated.


- I restarted with a new genome found here: https://www.internationalgenome.org/data-portal/sample/HG00114 . 
- Damn, downloading it sucked but I finally got it downloaded and expanded.
- Started with ERR013138_1.fastq
- Extracted 1000 records with seqkit.
- Extracted names (every 4th line after first) with awk (quality_scores.txt).
- Designed a regex to split out info from identifiers. *This regex is probably incorrect for some data formats.* : https://regex101.com/r/RnuEOJ/3

```bash
seqkit head -n 1000 ERR013138_1.fastq > thousand_records.fastq
awk 'NR % 4 == 1' thousand_records.fastq > identifiers.txt

``` 

- Parsed identifiers into a dictionary of component values (names, instrument, flow cell lane, etc). *Ibrahim: Does this look right?*

```bash
(bio) w134-87-145-178:flow_cell_lanes nat$ ./quap.py 
whole_identifier  : 
@ERR013138.1 IL39_4668:5:1:1035:1408/1
@ERR013138.2 IL39_4668:5:1:1035:8133/1
@ERR013138.3 IL39_4668:5:1:1035:6544/1
@ERR013138.4 IL39_4668:5:1:1035:2520/1
@ERR013138.5 IL39_4668:5:1:1036:15336/1

names  : 
@ERR013138.1
@ERR013138.2
@ERR013138.3
@ERR013138.4
@ERR013138.5

instrument  : 
IL39_4668
IL39_4668
IL39_4668
IL39_4668
IL39_4668

flow_cell_lane  : 
5
5
5
5
5

flow_cell_number  : 
1
1
1
1
1

x_coord  : 
1035
1035
1035
1035
1036

y_coord  : 
1408
8133
6544
2520
15336
```






## C++ Implementation

We need a working idea first, but eventually this will need implemented at production scale. The C++ source code for scalce is quite nice but
it has lots of extra features we might not need. I'm wondering what can we take from Scalce for our prototype?


## Not important- I tried to install Ibrahim's Benchmarking Toolkit 

Tried to follow the instructions https://github.com/sfu-compbio/compression-benchmark/blob/master/aws.md

, got as far as vagrant up before erroring out because of an "image not found" error.


# BLOCKAGES: Where we need help from Ibrahim

- Nat: I'm finding a hard time figuring out where exactly compression is happening in the Scalce source code (for reads).
- Nat: can Ibrahim compare and contrast, in Scalce, the compression method used for reads and the one used for quality scores? *Is our project essentially trying to adapt the loss-less method used for reads to quality scores, using flow-cell etc. info?*

- Nat: Sample data. I grabbed some human genome data from https://www.internationalgenome.org/data-portal/sample/HG00114 ; is this appropriate?



