# czid-prototype-ont-pipeline
Prototyping an ont-compatible metagenomic pipeline for czid. 

This is an active work-in-progress.

This repository contains the prototype .wdl implementation for initial testing of an ont-compatible metagenomics pipeline.

**Pipeline download and set-up:**

```
git clone git@github.com:katrinakalantar/czid-prototype-ont-pipeline.git
cd czid-prototype-ont-pipeline

docker build -t ontp ont-pipeline
```

**Running the pipeline on test data:**

Note, some of the set-up for this is currently done manually (i.e. building the minimap2 and centrifuge dbs). Detailed set-up notes (in rough form) are currently documented in `ont-pipeline/ont_pipeline_setup.txt`. These will be refined over time.


```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736-tiny.fastq minimap_host_db=reference/chr1.fa minimap_human_db=reference/chr1.fa NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent
```