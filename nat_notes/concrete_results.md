## Concrete Results

- Extract first n fastq records from a large file

```bash
seqkit head -n 1000 <large fastq file> > thousand_records.fastq
```

- Convert a large .fastq file to quality_scores.txt, identifiers.txt

```bash
awk 'NR % 4 == 0' <large fastq file> > quality_scores.txt
awk 'NR % 4 == 1' <large fastq file> > identifiers.txt
```

- Parse identifiers and make text files for analytics, output to a new directory identifierComponents/ . Overwrite everything in identifierComponents/ if it already exists.
- script is in scripts/ directory.

```bash
./parseidentifiers.sh identifiers.txt
```

- Convert quality scores to integers, write them to a text file, plot them

**This might be too slow for 3 GB full human genome files.**

```python
# In file quality_scores/quap.py
writeQualityScores("quality_scores.txt",'converted_quality_scores.txt')
plotFromConverted("converted_quality_scores.txt", 0, 100)
```

- Concatenating quality score data with flow cell lane, tile number, x coordinates, and y coordinates into a single text file, side_by_side.txt:

```bash
paste quality_scores.txt identifierComponents/flow_cell_lane.txt identifierComponents/flow_cell_number.txt identifierComponents/x_coord.txt identifierComponents/y_coord.txt > side_by_side.txt
```

- Extract average quality scores from a fastq file

With a python script from the scripts directory (might be too slow for large files):
```bash
python extract_average_qas.py <fastq file> <output file>
```

*Still fairly untested, but should work well. Let me know how it goes!*

With a bash script from the scripts directory (should work for all files including large ones):
```bash
chmod +x bash_extract_average_qas # Make script executable, needed the first time it's run
./bash_extract_average_qas.sh <fastq file> > <output file>
```

- Extract median quality scores from a fastq file

With a python script from the scripts directory (might be too slow for large files):
```bash
python extract_median_qas.py <fastq file> <output file>
```

*Still fairly untested, but should work well. Let me know how it goes!*
**Requires gnu awk, installed on macos with "brew install gawk"**
**It's a bit lazy about it's median calculation. Picks sorted_line[n/2] as the median even when n is even, when it really should be taking the average of sorted_line[floor(n/2)] and sorted_line[ceiling(n/2)]. Didn't seem worth it to do things better yet.**

With a bash script from the scripts directory (should work for all files including large ones):
```bash
chmod +x bash_extract_average_qas # Make script executable, needed the first time it's run
./bash_extract_lazy_median_qas.sh <fastq file> > <output file>
```

- Check read_numbers are simply increasing one by one, with c++:

This is much faster than the following bash method. From scripts/cpp directory, **making sure to hard code a filename in the source code**:

```bash
g++ -std=c++17 monotonic.cpp -out  # To build the executable
./out
```

With pure bash, from scripts directory:

```bash
chmod +x monotonic.sh
./monotonic.sh read_numbers.txt

```

- Get min, max values from a file with a single number on each line, **making sure to hard code a filename in the source code**:

With c++, from scripts/cpp directory:
```bash
g++ -std=c++17 minmax.cpp -o out
./out 

```





