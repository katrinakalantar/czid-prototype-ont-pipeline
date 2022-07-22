### RUNNING ACTUAL DATSETS THROUGH THE PIPELINE


**Launch a new, larger instance for building and running with full references**

https://czi.quip.com/ymxaAeuo0Fql/Individual-Attribution-Instances

```
./ec2ia launch --iam-role idseq-comp-bio --instance-type m5.4xlarge -- --availability-zone us-west-2c
```

Note, this instance is $0.75/hour = $18/day = $126/week

Use this instance to build host filtering databases.

Copy host genome files from s3 and make respective databases:

```
cd cd ONT-research/reference/

# HUMAN HOST - 
aws s3 cp s3://idseq-public-references/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_phiX_rRNA_mito_ERCC.fasta .
minimap2 -x splice -d hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi hg38_phiX_rRNA_mito_ERCC.fasta
minimap2 -x map-ont -d hg38_phiX_rRNA_mito_ERCC_mm-map-ont.mmi hg38_phiX_rRNA_mito_ERCC.fasta

# MOSQUITO HOST - 
aws s3 cp s3://idseq-public-references/host_filter/mosquitos/reference_fastas/mosquito_genomes_20181207.fa.gz .
minimap2 -x splice -d mosquito_genomes_20181207_mm-splice.mmi mosquito_genomes_20181207.fa.gz  # ERROR - segmentation fault
minimap2 -x map-ont -d mosquito_genomes_20181207_mm-map-ont.mmi mosquito_genomes_20181207.fa.gz  # ERROR - segmentation fault

# HUMAN REMOVAL - 
aws s3 cp s3://idseq-public-references/host_filter/human/2018-02-15-utc-1518652800-unixtime__2018-02-15-utc-1518652800-unixtime/hg38_pantro5.fa .
minimap2 -x splice -d hg38_pantro5_mm-splice.mmi hg38_pantro5.fa  
minimap2 -x map-ont -d hg38_pantro5_mm-map-ont.mmi hg38_pantro5.fa

```

Open the raw data from nanopore to get a .fastq file to use:
```
tar -xvf 20220615_ont_czid_fastqs.tar #there are 482 GB remaining space before unzipping the raw data.
```

Try to run host removal on a full sample on this instance. Bonus: see how far it gets into assembly!

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/1-idseq-hum.fq.gz minimap_host_db=reference/hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```


Stopped this at 27 hours into the assembly. It was running with --threads = 8 on a machine with 16 CPU. OH, it was running on an output of 1 million reads (instead of 100,000).  Prior to subsampling, there were 9066024 lines in the file = 2266506 sequences (2.2 million reads).


m5.16xlarge - 64 vCPU, so I'll update --threads to 64 in assembly step, $3/hour
```
./ec2ia launch --iam-role idseq-comp-bio --instance-type m5.16xlarge -- --availability-zone us-west-2c
```

Try to continue running assembly (on 100,000 reads)...
```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/1-idseq-hum.fq.gz minimap_host_db=reference/hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```
NOTE: results for 100k read subset of 1-idseq-hum is in: 20220720_232327_ontpipeline

Try to run 100k subsample of 2-idseq-hum
```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/2-idseq-hum.fq.gz minimap_host_db=reference/hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```
NOTE: results for 100k read subset of 2-idseq-hum is in: 20220721_181647_ontpipeline

Try to run 100k subsample of 3-idseq-hum
```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/3-idseq-hum.fq.gz minimap_host_db=reference/hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```
NOTE: results for 100k read subset of 3-idseq-hum is in: 20220721_211122_ontpipeline


aws s3 cp s3://idseq-prod-samples-us-west-2/comp-bio-workspace/ont-data/4-idseq-hum.fq.gz ../ONT-raw-data/
aws s3 cp s3://idseq-prod-samples-us-west-2/comp-bio-workspace/ont-data/5-idseq-hum.fq.gz ../ONT-raw-data/
aws s3 cp s3://idseq-prod-samples-us-west-2/comp-bio-workspace/ont-data/6-idseq-hum.fq.gz ../ONT-raw-data/


Try to run 100k subsample of 4-idseq-hum
```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/4-idseq-hum.fq.gz minimap_host_db=reference/hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```
NOTE: results for 100k read subset of 4-idseq-hum is in: 20220721_234114_ontpipeline

Try to run 100k subsample of 5-idseq-hum
```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/5-idseq-hum.fq.gz minimap_host_db=reference/hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```
NOTE: results for 100k read subset of 5-idseq-hum is in: 20220722_022439_ontpipeline


Try to run 100k subsample of 6-idseq-hum
```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/6-idseq-hum.fq.gz minimap_host_db=reference/hg38_phiX_rRNA_mito_ERCC_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```
NOTE: results for 100k read subset of 6-idseq-hum is in: 20220722_051453_ontpipeline



TRYING A MOSUITO SAMPLE EVEN THOUGH THERE WAS A SEG FAULT IN GENERATING THOSE HOST GENOMES -- 
aws s3 cp s3://idseq-prod-samples-us-west-2/comp-bio-workspace/ont-data/17-idseq-mosq.fq.gz ../ONT-raw-data/

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=../ONT-raw-data/17-idseq-mosq.fq.gz minimap_host_db=reference/mosquito_genomes_20181207_mm-splice.mmi minimap_human_db=reference/hg38_pantro5_mm-splice.mmi library_type=RNA NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```
WAS RUNNING ALIGNMENT BUT HAD ERROR FOR --split-prefix, so I stopped it. The issue right now is space on the instance. 


I may need to have a function to copy files into an output directory so that I can remove the unnecessary files (?) 


I may need to zip the results from the pipeline runs and move them elsewhere to free up space for running new pipelines!! 
i.e. s3://idseq-prod-samples-us-west-2/comp-bio-workspace/

tar -czvf 20220722_051453_ontpipeline.tar.gz 20220722_051453_ontpipeline
aws s3 cp 20220722_051453_ontpipeline.tar.gz s3://idseq-prod-samples-us-west-2/comp-bio-workspace/
rm 20220722_051453_ontpipeline.tar.gz
rm -rf 20220722_051453_ontpipeline
