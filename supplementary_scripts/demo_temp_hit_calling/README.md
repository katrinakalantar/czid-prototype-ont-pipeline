This directory contains prototype of hit-calling script that can be used to generate interim taxonomic relative abundances (in reads and bases at genus and species level).

To run:
```
python3 tally_hits_all.py 20 20220823_145511_ontpipeline.reads-to-contigs.out
```

Note: the `20220823_145511_ontpipeline.reads-to-contigs.out` file is a subset of the reads-to-contigs output file from the pipeline (columns 1, 3, 10).

Note: 20220823_145511_ontpipeline.reads-to-contigs.out contains just the top 1000 lines of the original file due to file size limitations