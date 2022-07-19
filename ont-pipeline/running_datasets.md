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