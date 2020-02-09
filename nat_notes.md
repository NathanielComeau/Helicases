# Nat's Semi-Private Notes

Just some ill-formed notes to myself

## Nice to haves?

- Convert a large .fastq file to quality_scores.txt, identifiers.txt

```bash
awk 'NR % 4 == 0' <large fastq file> > quality_scores.txt
awk 'NR % 4 == 1' <large fastq file> > identifiers.txt
```

- Convert quality scores to integers, write them to a text file, plot them

```python
# In file quality_scores/quap.py
writeQualityScores("quality_scores.txt",'converted_quality_scores.txt')
plotFromConverted("converted_quality_scores.txt", 0, 100)
```
- Parse identifiers and make text files for analytics

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









