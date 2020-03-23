# Lab Notebook

This is my notebook where I write down what I tried so other people can follow in my footsteps.

## Playing with scalce and initial explorations

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

### Attempts to extract identifiers and qa from full human genome

- Tried to do it using Python; failed when it tried to use 12GB of RAM. Although it wrote 600,000 of the identifiers before failing; maybe some small optimizations are all that's needed to do it straight from Python?

```python
# In file quality_scores/quap.py
writeQualityScores("quality_scores.txt",'converted_quality_scores.txt')
plotFromConverted("converted_quality_scores.txt", 0, 100)
```

```python
```python
# Setup
import sys
import re
import os
import shutil

def splitIdentifiers(idFile):
    with open(idFile, 'r') as f:
        identifiers = f.readlines()
    identifiers = [x.rstrip('\n') for x in identifiers]

    # example id: @ERR013138.1 IL39_4668:5:1:1035:1408/1
    # Splits into name, instrument, flow_cell_lane, flow_cell_number, x_coord, y_coord
    # Regex pattern1 based on https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers: 
    #               (\@\w*.\d+)\s(\w+):(\d+):(\d+):(\d+):(\d+)(#\d+)*\/(\d+)$
    #               https://regex101.com/r/RnuEOJ/3
    p = re.compile(r'(\@\w*.\d+)\s(\w+):(\d+):(\d+):(\d+):(\d+)(#\d+)*\/(\d+)$')
    id_components = {"whole_identifier":[],
              "names":[], 
              "instrument":[], 
              "flow_cell_lane": [],
              "flow_cell_number": [],
              "x_coord": [],
              "y_coord": [],
              #"index_number": [],
              "pair_member": []
              }

    for ident in identifiers:
        #print(ident)
        match = p.match(ident)
        if not match:
            print("ERROR: Pattern matching failed on identifiers!")
            sys.exit()
        id_components["whole_identifier"].append(match.group(0))
        id_components["names"].append(match.group(1))
        id_components["instrument"].append(match.group(2))
        id_components["flow_cell_lane"].append(match.group(3))
        id_components["flow_cell_number"].append(match.group(4))
        id_components["x_coord"].append(match.group(5))
        id_components["y_coord"].append(match.group(6))
        #id_components["index_number"].append(match.group(7))
        id_components["pair_member"].append(match.group(8))
    
    # print first 5 id components to make sure we did it right
    # for key in id_components:
    #     print(key, ' : ')
    #     for i in range(5):
    #         print(id_components[key][i])
    #     print()
    return id_components

def writeIdentifiers(id_components):
    # Write components of identifiers to individual text files in a new directory 
    # Make a new directory to hold our files
    dirName = 'identifierComponents'
    try:
        os.mkdir(dirName)
    except FileExistsError:
        shutil.rmtree(dirName)
        os.mkdir(dirName)

    for key in id_components:
        with open('identifierComponents/'+key+'.txt', 'w') as f:
            for item in id_components[key]:
                f.write(item+'\n')

def reassembleIdentifiers():
    # Just a quick sanity check that I parsed identifiers properly
    id_components = {"whole_identifier":[],
              "names":[], 
              "instrument":[], 
              "flow_cell_lane": [],
              "flow_cell_number": [],
              "x_coord": [],
              "y_coord": [],
              #"index_number": [],
              "pair_member": []
              }
    for key in id_components:
        with open('identifierComponents/'+key+'.txt', 'r') as f:
            lines = f.readlines()

        lines = [x.rstrip('\n') for x in lines]
        id_components[key] = lines

    with open('reassembled.txt', 'w') as f:
        for i in range(len(id_components["names"])):
            # example id: @ERR013138.1 IL39_4668:5:1:1035:1408/1
            f.write(f'{id_components["names"][i]} {id_components["instrument"][i]}:{id_components["flow_cell_lane"][i]}:{id_components["flow_cell_number"][i]}:{id_components["x_coord"][i]}:{id_components["y_coord"][i]}/{id_components["pair_member"][i]}\n')

result = splitIdentifiers('identifiers.txt')
writeIdentifiers(result)
# reassembleIdentifiers() # Can do this to check parsing didn't mess up the files
```
```

- Maybe we can do it in bash somehow?

- Fuck me, that was annoying. Doing it in awk seems to use a constant amount of memory (~ 400 KB):

```bash
awk -F'[ .:#/]' '{print $1, $2, $3, $4, $5, $6, $7, $8}' identifiers.txt

# This doesn't work for some reason...?
#awk -F'[ .:#/]' '{print $1 > "read_name.txt" print $2 > read_number.txt}' identifiers.txt

# Fine, we'll do it iteratively the old fashioned way:
echo "Starting names"
awk -F'[ .:#/]' '{print $1}' identifiers.txt > read_names.txt
echo "Starting numbers"
awk -F'[ .:#/]' '{print $2}' identifiers.txt > read_numbers.txt
echo "Starting instruments"
awk -F'[ .:#/]' '{print $3}' identifiers.txt > instrument_names.txt
echo "Starting lanes"
awk -F'[ .:#/]' '{print $4}' identifiers.txt > flow_cell_lanes.txt
echo "Starting tile numbers"
awk -F'[ .:#/]' '{print $5}' identifiers.txt > tile_numbers.txt
echo "Starting x_coords"
awk -F'[ .:#/]' '{print $6}' identifiers.txt > x_coords.txt
echo "Starting y_coords"
awk -F'[ .:#/]' '{print $7}' identifiers.txt > y_coords.txt
# We could try to grab the index number / paired end info as well, fuck that for now 
#awk -F'[ .:#/]' '{print $8}' identifiers.txt > paired_ends.txt

```
- Testing it on whole genomes. Started at 11:56 am. Only using 388 kb of memory.
- Fuck yeah!! Took about a minute to finish one field!! We're in business. Took ~ 10, 15 minutes to finish.
- Gonna polish it into a shell script that spits out fields into a new directory, identifierComponents, and is used like:

```bash
./parseidentifiers identifiers.txt
```

```bash
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
```


### Converting large amounts of quality scores to integers

- Not sure how. Did some bash/sed googling and found nothing easy. Do we even need this?










### Analysis of Preliminary Identifier Files

- Extracted identifier fields from "ERR013138_1.fastq" with awk
- 30,003,628 identifiers in total

```bash
In [6]: ls -lh                                                                  
total 2684616
-rw-r--r--  1 nat  staff    57M 12 Feb 12:37 flow_cell_lanes.txt
-rw-r--r--  1 nat  staff   286M 12 Feb 12:35 instrument_names.txt
-rw-r--r--  1 nat  staff   315M 12 Feb 12:32 read_names.txt
-rw-r--r--  1 nat  staff   247M 12 Feb 12:34 read_numbers.txt
-rw-r--r--  1 nat  staff    89M 12 Feb 12:39 tile_numbers.txt
-rw-r--r--  1 nat  staff   158M 12 Feb 12:41 x_coords.txt
-rw-r--r--  1 nat  staff   159M 12 Feb 12:42 y_coords.txt
```

- flow_cell_lanes.txt: everything equaled 5.
- instrument_names.txt: everything equaled 'IL39_4668'.
- read_names.txt: everything equaled '@ERR013138'
- read_numbers.txt: just a bunch of indices from 1 to the total number of reads (30,003,628)

- *tile_numbers.txt*: Holy shit, it's a quadratic or something!! See tile_number_frequencies1.png or tile_number_frequencies2.png for reference.

- x_coords: **Finish this**
- y_coords: **Finish this**

- Both x and y plotted together:
    - HA! First try, it turned into a huge brick. Meaning I think they look to be pretty uniformly distributed over a rectangular range.
    - Second attempt showed distribution to be a rounded edge rectangle, was still too slow to really manipulate it.
    - Let's see if it's uniform over the x and y range? *However, I'm not entirely convinced uniformity in the frequency diagrams means uniformity in the distributions of the data.* Let's leave this for now.

    - HEYYO, let's try to do it with Gaussian filtering!!

    - Got some pretty results:

- **IDEA: Plot average quality score as a function of x and y!!!!!**

### Plotting quality scores as a function of x and y

- Yo, maybe let's start with just 1000 reads to keep things simple.

```bash
seqkit head -n 100000 ERR013138_1.fastq > little.fastq

```
_ *Huh, turns out read lengths are not always 100 bp? I was assuming that, weird*

- I spot checked some of the quality scores by plotting them with qa_manipulation.py. They maybe had a little bit of locality when I plotted 10 at a time randomly sampled. *Maybe there's potential for ordering similar quality scores with a little bit of loss to make them highly compressible?* Eg, take all quality scores that slope downward to the right in a similar way, group them (fuzz with a little bit of loss?) and try to compress them.

- Calculated average quality scores and saved them in a text file as 'average_quality_scores.txt'

- Extracted identifier fields and moved x and y coordinates into appropriate directory

```bash
cp ../../scripts/parseidentifiers.sh ./
./parseidentifiers.sh identifiers.txt 
mv identifierComponents/x_coords.txt ./
mv identifierComponents/y_coords.txt ./
```

- Holy shit, not sure where this is going. See x_and_y_explorations.ipynb for details.

- **Need to come back and actually plot quality scores as a function of x and y**

- Maybe I can weight x,y points by multiplying them by their average quality score? So if a point has a quality score of 16, we'd add 16 copies of that score to the diagram...

- Tried to do this in x_y_and_average_qa_exploration.ipynb. Didn't really see much in the plots. I think differences in x and y density are wiping out differences in quality scores. Trying to normalize x and y in two dimensions feels really messy. Picture a 3d field of points (x, y, quality score) as well. I feel like with all the smoothing we're doing (averaging, binning, gaussian smoothing), we really won't pick up on correlations of qa and x, y before doing some clustering.

- **Try 3d surface map, topographical map?**


- Made a figure trying to subtract one from the other; gotta interpret that at some point. **Try doing it again for the whole human genome**

- [Binned statistics were pretty cool! Both 1D and 2D](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binned_statistic.html) . With ERR013138_1.fastq, I extracted x, y, and quality scores and used the scipy binned statistics modules to create some figures exploring average quality score as a function of x and y. See big_x_y_qa_analysis.ipynb for details. Made 3 figures: 

    - full_genome_x_and_y_and_qa_average_binning_smoothing.png
    - full_genome_x_and_y_and_qa_medianing_binning_smoothing.png
    - full_genome_x_and_y_counts_binning.png

- I binned the 100,000 records file, wasn't really looking good for 1000 bins but made a figure with various sigma levels and 10,50,100 bins. Saved the figure in figures. Gaussian smoothing seems to be doing weird shit to it for small numbers of bins.


### Analysis of Ibrahim's large (10-17GB) human genome fastq files

- Downloaded ERR174325_1.fastq.gz from [deez](https://github.com/sfu-compbio/compression-benchmark/blob/master/samples.md)
    - HA! I didn't even have enough disk space on my laptop to decompress it. Let's get access to a new machine soon...

- Trying with a slightly smaller file (11GB) ERR174340_1.fastq.gz, hopefully this will work.

- Decompressed with gzip -d ERR174340_1.fastq.gz -c > ERR174340_1.fastq
    - **Failed after unpacking 19 GB.**
    - Tried decompressing to my hard drive; it worked, but only after unpacking 36 GB!!

- Extracted quality scores and identifiers
    - started quality scores at 10:38am: awk 'NR % 4 == 0' /Volumes/Marvin/ERR174340_1.fastq > quality_scores.txt
    - Finished at 11:16am, with final size of 14GB
    - started identifiers at 11:19am with awk 'NR % 4 == 1' /Volumes/Marvin/ERR174340_1.fastq > identifiers.txt.
    - Finished at 11:56am, with final size of 7.3GB


- Parsed identifiers with "./parseidentifiers.sh identifiers.txt"
    - Started around 1:15pm, took until 2:30pm to finish

- Analysis of identifier fields
    - **Did I even split these apart correctly? Need to double check that my assumptions for the  regex are correct (that the regex is semantically correct).**
    - Read names seemed to all be the same... Leaving those aside for now.
    - Same for instrument names... Leaving those aside for now.
    - Read numbers are monotonically increasing; I checked with a cpp and also with a bash script. Leaving them be for now.
    - flow_cell_lanes: every single value is 3. Leaving them for now.
    - tile_numbers: seems to be monotonically increasing as a function of read_number. [Here's the output from searching linearly for the max from the start to the end of tile_numbers.txt](snips/freq.md). There might be interesting stuff in here; I should plot frequencies and tile_number as functions of read_number, x_coord, and y_coord at some point.
        - The max was 2308, min 1101.
        - Frequencies looked like [this, raw.](snips/2.md)
        - [As an image.](figures/tile_number_frequencies.png)

    - x_coords: max was 21403, min was 1000.
        - Frequencies looked like [this, raw.](snips/3.md)
        - [As an image.](figures/x_coord_frequencies.png)
    - y_coords: max was 200940, min was 2073.
        - Frequencies looked like [this, raw.](snips/4.md)
        - [As an image.](figures/y_coord_frequencies.png)


- Extraction of quality scores
    - File is 14GB in size. Lettttt's try it with bash_extract_average_qas.sh ?? Might not be worth it to rewrite in c++...
        - Would have taken 9 hours with the shell script!! Time to try in C++
        - **Huh, nevermind, it took 40 minutes. Hahahahaha**
        - **Stopping here because I'm sick of data wrangling for now**



- Mosaicing x and y coordinates

    - For plotting quality score as a function of x and y, we know my old system works. It worked to plot ~ 30,000,000 x, y, and qas together.
    - This new data has ~ 150,000,000 datapoints (x, y, and qa). It's only ~ 5 times larger!! So if we bin the data into 5 bins, we should be good. 
    - Based on how the frequency plots look, lets bin in y this time.
    - Range is 198,867, into 5 bins: 39,773.4 values per bin. So the ranges will be 
        - [2072, 41846, 81619, 121392, 161165, 200941]
    - Need to sort x, y, qa by the y coordinate.

        - x, y coord files are 1GB in size.
        - BUT the resulting bins don't have to be sorted. We just need to bin the data.
        - One pass through both x and y file, if y meets condition, write the point out to the relevant bin file.
        - Requires keeping 3 + number_bins files open at same time.
        - Wrote a barebones c++ script to do it, used 380kb of memory (!!!!), took ~ 20 minutes to run.
            - Checked that the [file sizes added up.](snips/5.md)
            - **Not super sure how to check it's right. Might have fucked this up here.**
            - **Stopping to get quality scores, and do it again with quality scores.**














*Maybe I'll have to mosaic to be able to make the plots?* That would make things much, much easier than trying to devise a clever way to load things into memory all at once.

### K-means clustering of quality scores

So far I have 3d data (x, y, average quality score). Nothing super obvious in terms of corelations popped up when I tried to plot it but I might have been binning super shittily (probably?). 

- Started exploration with first 100,000 x-coords, y-coords, average qas.



















# Nat's Notes to Himself

Just some ill-formed notes to myself


## Big Picture

I think I'm starting to see what the big picture is.
    - Use a different model for a particular x,y/flow cell lane/tile "chunk". If we can cluster quality scores by those chunks effectively, we're golden.

To figure out how to chunk, I think we need to do the binning analysis again on some of the test datasets Ibrahim provided.

Following that (or possibly in tandem) we'll need to try clustering this data with something like k-means. **Question: how else can we try and cluster this data?**

Efficient implementation will be a pain (?).

Paper writing will be cool. Good to see there is similar literature out there.



## Nice to Haves




## Immediate ToDos

- **Is my parsing (my regex) of identifier fields correct?**

- Redo binning and plotting of average quality score as a function of x and y for Ibrahim's data.

- Get a basic compressor up and running. Maybe we could try to reproduce Ibrahim's results for his simplest model? Eventually we'll want to start running that thing with different models we feed to it.

- Mayyyybbbeee try to visualize data as a data cube / butterfly plot?

- Try to understand arithmetic coders a little better, including (?) the one used by Samcomp. I would like to build a simple one that reads two streams, data and 1's or 0's, that switches between two models based on the 1's and 0's stream.




- Do some more fun plotting of x and y coordinates with filters. Try to figure out if x and y are truly uniform (possibly a 2d standard deviation or something?). Nearest neighbours?

- Finish exploring x_coords.txt and y_coords.txt

- Essentially we want to see how those quality scores correlate with identifier fields.


### Shelved

- Try to plot x and y, with colour representing average quality score. Do it again for median quality score. Maybe you'll have to try a weird binning/filtering thing with it.

	- This felt really messy; I might need to find a way to normalize in x and y to do this right. There probably is a good way to visualize this 3d data (x, y, average quality score), I'm just not quite sure how to bin it yet.

- Got some good identifier field extractions in /Users/nat/research/Helicases/scratch/identifierComponents . Why not run some basic analytics on each of those text files? See what values they take on, frequencies, averages, periodicities, etc.

    - **Why is the frequency count distorted for the tile numbers?** See tile_number_frequencies1.png or tile_number_frequencies2.png
    - Are the x and y coordinates truly uniformly distributed in a rounded rectangle? I could make a big plot (takes time, doesn't find subtle non-uniformity). I could check that the average distance between points equals the expected distance between points (fucking hard to do computationally, might not catch small variations either). Maybe make a big-ass plot and go from there... How to figure out of something is uniformly distributed?



## Eventual things to get around to

- Look up papers on quality score recalibration, maybe something in there is useful to us.

- What if we plot flow cell lane/ tiles with x and y? Can that explain some of the structure we see in quality score as a function of x and y?

- k-means clustering
- Which identifier components should we use? Which can we throw out, they're useless?
    - Probably can throw out cell lanes, instrument names, read names, most of read numbers.
    - Only interesting values are x, y, flow tile. 
- Do a quick sanity check on identifier components to explore their properties.
    1) Is name & instrument the same for all identifiers?
    2) Are all their lengths equal?
- Understand move to front, run length, arithmetic encoding.


## Done

- Try to plot the full human genome x, y, quality score
    - Done, need to redo it for Ibrahim's samples. Does it scale to 17 GB?


## Crazy Ideas and Notes



- Try binning by tile number. Looks like reads are roughly distributed uniformly between tiles. Maybe qas are correlated within each tile.


- **Compression Ideas**
  - There are vertical streaks in the quality score x y images. What if each x-coordinate has a different model?
  - Highly binned data has differing quality scores. What if each bin had it's own model?
  - *I really don't know how tile numbers fit into this. We shouldn't forget about them.*


- How about we stack multiple x,y,qa images together to try and detect systematic biases in the illumina platform?

- Turn each set of x, y, quality scores into a data cube!!! Where the 3rd axis is quality score.
    - Are x and y coordinates unique? I'm just wondering how I'll convert everything to a rectangular grid to make a data cube out of them.

- What does the best possible input to gzip look like? Can we find an optimal ordering of reads, 

- Can we use NA12878 and GIAB to validate any of these results?


- Crazy data vizualization ideas:

  - [Cool contour plots](https://plot.ly/python/3d-surface-plots/)
  - Try to do some kind of surface plots?
  - Voronoi diagrams

  - [Awesome Gaussian filtering to determine nearest neighbors](https://stackoverflow.com/questions/2369492/generate-a-heatmap-in-matplotlib-using-a-scatter-data-set)
  - [More general issues plotting large data sets](https://stackoverflow.com/questions/4082298/scatter-plot-with-a-huge-amount-of-data)


- K-means clustering ideas
    - [How to pick k for k-means](https://www.quora.com/How-can-we-choose-a-good-K-for-K-means-clustering)
    - [scikit-learn k-means](https://scikit-learn.org/stable/modules/clustering.html#k-means)
    - [Voronoi diagrams?](https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.spatial.Voronoi.html)



- Cool plotting of mosaicing figures https://python4astronomers.github.io/plotting/advanced.html











