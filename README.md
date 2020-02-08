# Helicases

# Nat's work on weekend of Feb 8th

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



