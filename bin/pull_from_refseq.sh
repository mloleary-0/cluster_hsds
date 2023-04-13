#!/usr/bin/env bash
DOWNLOAD_DEST_DIR=(refseq_genomic_fnas)
GENBANK_TABLE_FILE=$(realpath config/prokaryotes_11042022.csv) # add this to $MY_RUN_DIR


mkdir -p $DOWNLOAD_DEST_DIR
cd $DOWNLOAD_DEST_DIR
while IFS="," read OrganismName OrganismGroups Strain BioSample BioProject Assembly Level SizeMb GC Replicons WGS Scaffolds CDS ReleaseDate GenBankFTP RefSeqFTP etc
do
GenBankFTPadjust=$(echo $RefSeqFTP | sed 's/"//g')
TARGET=($GenBankFTPadjust'/'*'_genomic.fna.gz')
echo $TARGET
wget --reject "GCF_006369915.1_ASM636991v1*","GCF_006369965.1_ASM636996v1*","GCF_000166855.2_ASM16685v2*","GCF_000767565.1_ASM76756v1*","GCF_006370085.1_ASM637008v1*","GCF_006370035.1_ASM637003v1*","GCF_006370055.1_ASM637005v1*","GCF_006369895.1_ASM636989v1*","GCF_006369905.1_ASM636990v1*","GCF_006369995.1_ASM636999v1*","GCF_006370045.1_ASM637004v1*","GCF_001572105.1_ASM157210v1*","*_from_*" $TARGET
done < $GENBANK_TABLE_FILE
cd ..


#--- get x. taiwanensis genomes, refseq version.
# PLS299: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/177/435/GCF_013177435.1_ASM1317743v1
# PLS235: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/352/785/GCF_003352785.1_ASM335278v1
# PLS244: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/910/935/GCF_009910935.1_ASM991093v1
# note there is a bad PLS229 in here, which is the fourth sequenced strain.
cd $DOWNLOAD_DEST_DIR
for FTP in ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/177/435/GCF_013177435.1_ASM1317743v1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/352/785/GCF_003352785.1_ASM335278v1 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/910/935/GCF_009910935.1_ASM991093v1 
do
wget --reject "*_from_*" $FTP'/'*'genomic.fna.gz'
done
cd ..

# # for all x. taiwanensis genomes, use:
# cd $DOWNLOAD_DEST_DIR
# while IFS="," read OrganismName OrganismGroups Strain BioSample BioProject Assembly Level SizeMb GC Replicons WGS Scaffolds CDS ReleaseDate GenBankFTP RefSeqFTP etc
# do
# GenBankFTPadjust=$(echo $RefSeqFTP | sed 's/"//g')
# TARGET=($GenBankFTPadjust'/'*'_genomic.fna.gz')
# echo $TARGET
# wget --reject "GCF_000576405.1","*_from_*" $TARGET
# done < $GENBANK_TABLE_FILE
# cd ..



#------------ Get Dixon specifically, since it is suppressed:
cd $DOWNLOAD_DEST_DIR
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/166/835/GCF_000166835.1_ASM16683v1/GCF_000166835.1_ASM16683v1_genomic.fna.gz # Dixon, since it is suppressed
cd ..
