# cluster_hsds
A snakemake pipeline to extract Type I RM system components and cluster hsdS genes from Xylella fastidiosa.  Simply, this will check all nucleotide CDS prediction files with the `.fna` extension in a folder `genomes` for all Hsd components belonging to XfaI-IV (e.g., HsdM, HsdR, HsdS), cluster HsdS alleles with type allele references using `cd-hit`, and use `make_multi_seq.pl` to split these clusters out into seperate fasta-formatted files.  See dag.pdf for the specific workflow.
\
# Setup 
## Environment Setup 
Pull down the repository into your run directory using:
```
https://github.com/mloleary-0/cluster_hsds.git
``` 

Set up the conda environment with `ncbi-blast`, `cd-hit`, and `snakemake` using the provided `environment.yml` file: \
```
conda env create -f config/environment.yml
``` 
or 

```
mamba env create -f config/environment.yml
``` 

## Data Setup 
If you have a set of fasta-formatted nucleotide CDS ready to go, make a folder called `genomes` and move your genomes into it.  Note that by default, this pipeline acts on files with the `.fna` extension, but that can be overridden to use another extension (see Options below).  I chose `.fna` as default because this is the format the relevent GenBank and RefSeq data files come in. 

```
mkdir genomes
mv sample_genome.fna genomes/
```

Alternatively, you may want to fetch CDS files from RefSeq directly.  For example, you could fetch all Xylella fastidiosa genomes from RefSeq to `genomes` using the NCBI `prokaryotes.txt` file.  This file is downloadable from the [Genome page / Genome Assembly and Annotation report page](https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/173/).  Unfortunately I am not sure where a direct link to it is is, but you can manually download it by clicking on the "Download" icon above the table.  Then, you could run something like the code block below to pull all assemblies with a RefSeq entry down into `genomes` and unzip them: 

```
A=(genomes) 
B=$(realpath prokaryotes.csv) 
mkdir -p $A && cd $A 
while IFS="," read OrganismName OrganismGroups Strain BioSample BioProject Assembly Level SizeMb GC Replicons WGS Scaffolds CDS ReleaseDate GenBankFTP RefSeqFTP etc ; do
  GenBankFTPadjust=$(echo $RefSeqFTP | sed 's/"//g') 
  TARGET=($GenBankFTPadjust'/'*'_genomic.fna.gz') 
  echo $TARGET 
  wget $TARGET 
done < $B 
cd .. 
gunzip $A/*.gz
```
If you wanted to use the same set of RefSeq genomes I used in my 2022 manuscript for some reason, run:
```
chmod +x bin/pull_from_refseq.sh
bash bin/pull_from_refseq.sh
```
Note that in this case, I rejected some strains that were duplicates, labortory mutant strains, or appeared to have other issues. \

# Running the pipeline 
Now that you have your CDS nucleotide predictions in the `genomes` folder in `.fna` format, running the pipeline is very easy.  Activate the conda environment, and launch the pipeline with `N` threads: 
```
conda activate hsd_cluster_env
snakemake -jN
```
## Options
Note that blast e-value, the input directory, and .fna file extension can be overwritten in config/pipeline_config.yaml, or using the command line options shown below (default values given for example):
```
snakemake -j N --config evalue=1e25 extension=.fna genomes_dir=genomes
```

# Output
An overview of the intermediate steps and output: 
```
01-blastdbs 
02-blastn 
03-blastdbcmd 
05-CD-HIT 
  inputs
  output
    XfaI
      XfaI.clusters # <- clusters file with representative sequences for each cluster 
      XfaI.clusters.clstr #<- an overview of the clusters 
      multi-seq
        0 #<- fasta-formatted file of every sequence in cluster 0
        1 #<- fasta-formatted file of every sequence in cluster 1
        2 #<- fasta-formatted file of every sequence in cluster 2
        ...
     ...
 ```
 ## Quality control
Note that this is not foolproof - in my experience, truncated HsdS genes will end up in mixed clusters due to the conserved regions, or sharing a TRD but not the other.  If you see something that is shorter than expected (say, under 1.1 kb), I strongly recommend doing a manual curration.  One disadvantage to the method I'm using here is that I rely on a gene prediction rather than going to the chromosome sequence itself, where I would be able to extract the entire hsdS gene even if it has a frameshift or other mutation.  However, it is difficult to identify HsdS alleles with novel TRDs using that approach, which is why I rely on CDS here.
