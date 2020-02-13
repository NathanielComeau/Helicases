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

**This is too slow for 3 GB full human genome files. Only useful for small samples.**

```python
# In file quality_scores/quap.py
writeQualityScores("quality_scores.txt",'converted_quality_scores.txt')
plotFromConverted("converted_quality_scores.txt", 0, 100)
```

- Concatenating quality score data with flow cell lane, tile number, x coordinates, and y coordinates into a single text file, side_by_side.txt:

```bash
paste quality_scores.txt identifierComponents/flow_cell_lane.txt identifierComponents/flow_cell_number.txt identifierComponents/x_coord.txt identifierComponents/y_coord.txt > side_by_side.txt
```


