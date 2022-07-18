
# Developing Initial ONT PIPELINE

### Steps for incremental development:

Build the docker image and run very first pipeline run [bare-bones]
```
docker build -t ontp ont-pipeline/
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq 
```
Set up the host genome
```
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
```

note: not creating an actual index now, but instead just calling minimap2 on the .fa file
GRCh38_latest_genomic.fna 

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=GRCh38_latest_genomic.fna
```

 **oops**, this is running out of memory on my small test instance. So, I just downloaded CHR1 for testing.
```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=reference/chr1.fa
```

**Added human filtering and subsampling steps!**

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=reference/chr1.fa minimap_human_db=reference/chr1.fa
```

**Adding flye for assembly**

If you run into a space image during docker build with error message `You don't have enough free space in /var/cache/apt/archives/`, then run the following:
```
docker container prune; docker image prune; docker volume prune
```

**adding in alignment step**
Now have two additional parameters to specify the databases.

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=reference/chr1.fa minimap_human_db=reference/chr1.fa NT_minimap2=reference/dummy_nt_minimap2 NT_centrifuge=reference/dummy_nt_centrifuge alignment_test_mode=all_mm
```

**Downloading the refseq database** for bacteria,viral,fungi
Using the `ncbi-genome-download` tool available here: https://github.com/kblin/ncbi-genome-download
```
ncbi-genome-download --refseq-categories reference --formats fasta bacteria,viral,fungi   
# ^^ this was taking TOO LONG
ncbi-genome-download --assembly-levels complete --parallel 4 --formats fasta bacteria,viral,fungi
# ^^ still too slow / clunky for this rapid testing prototype
```

Just using files I had previously retrieved (though not sure where I got them)...
```
cat ../reference_files/bacteria.1.1.genomic.fna.gz ../reference_files/viral.1.1.genomic.fna.gz > bacterial_viral_dummy_db.fna.gz
```

**Build NT database for minimap2**
```
minimap2 -x asm20 -d mm-asm20_bacterial_viral_dummy_db.mmi bacterial_viral_dummy_db.fna.gz 
```

Now, running the workflow with the actual database
```
using actual minimap2 index
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=reference/chr1.fa minimap_human_db=reference/chr1.fa NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/dummy_nt_centrifuge alignment_test_mode=all_mm
```



**Adding Centrifuge to the NT alignment step**
Requires update to the Dockerfile

https://github.com/DaehwanKimLab/centrifuge/archive/refs/tags/v1.0.4.tar.gz

Downloading centrifuge in a way that can be replicated in Dockerfile:
```
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip
unzip centrifuge-1.0.3-beta-Linux_x86_64.zip
```

Downloading the reference files (manually)
```
cd ./tools/centrifuge-1.0.3-beta/
./centrifuge-download -o taxonomy taxonomy
centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map 
# ^^ this command didn't work, so following instructions here: https://github.com/DaehwanKimLab/centrifuge
cd ./centrifuge-1.0.3-beta/indices
make p+h+v
# requires installation of BLAST+
sudo apt install ncbi-blast+
make p+h+v # still fails
```

I'm having trouble creating a new centrifuge index using the linux binary version, so trying to use the pre-existing index that I generated using the version of centrifuge that I compiled locally from github (instead of the linux binary version). The index is here: ~/comp-bio-secure/ONT-research/tools/centrifuge/indices/reference-sequences/p_compressed+h+v

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=reference/chr1.fa minimap_human_db=reference/chr1.fa NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=tools/centrifuge/indices/reference-sequences/p_compressed+h+v alignment_test_mode=split_mm_cent
```

File not found, need to create a "file" of the database directory and then unzip that for use in the alignment. Re-zipped the directory and updated the wdl to reflect this.

```
cp -r tools/centrifuge/indices/reference-sequences/ ./reference/centrifuge-ref
zip -r reference/centrifuge-ref.zip reference/centrifuge-ref
```

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=reference/chr1.fa minimap_human_db=reference/chr1.fa NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent
```

**Making DIAMOND database**

DIAMOND / NR alignment requires protein sequences. Downloading the NR database.
```
cd reference
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz # stopped this short to get just a tiny .fa file (1% of the db)
mv nr.gz tiny-nr.gz

# continued this command in screen to eventually get the full db [detached from 25807.pts-2.ia-katrina-kalantar-1656703242]

# building the DB
./tools/diamond makedb --in ./reference/tiny-nr.gz -d ./reference/tiny-nr
```

Running the wdl pipeline with diamond db parameter

```
miniwdl run --verbose ont-pipeline/run.wdl docker_image_id=ontp input_fastq=data/SRR12458736.fastq minimap_host_db=reference/chr1.fa minimap_human_db=reference/chr1.fa NT_minimap2=reference/mm-asm20_bacterial_viral_dummy_db.mmi NT_centrifuge=reference/centrifuge-ref.zip alignment_test_mode=split_mm_cent NR_diamond=reference/tiny-nr.dmnd
```


**Adding medaka capability / options**


We want to get the medaka option running! It failed to build the db based on docker space issue.

Update the docker location to enable building with medaka - using instructions here: https://github.com/katrinakalantar/vigilant-palm-tree/blob/main/2021_0315/notes_week.md#how-to-update-the-size-of-the-storage-on-the-instance-so-docker-will-run

Try to build docker image with medaka line uncommented.





