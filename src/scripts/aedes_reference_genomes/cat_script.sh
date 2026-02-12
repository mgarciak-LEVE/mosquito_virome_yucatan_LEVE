#!/bin/bash
eval "$(/Users/Parsimony/miniconda3/bin/conda shell.bash hook)"

echo "=== Building Aedes Superreference  for Hodt Depletion ==="

# Script for Aedes genome and chromosome annotation concatenation
# It is recommended to prefix chromosome names. This way we avoid naming conflicts

# Aedes aegypti (RefSeq, chromosome-level)
sed 's/^>/>Aaeg_|/g' GCF_002204515.2_AaegL5.0_genomic.fna > Aedes_aegypti.fna
# Aedes albopictus (RefSeq, chromosome-level, 2024)  
sed 's/^>/>Aalb_|/g' GCF_035046485.1_AalbF5_genomic.fna > Aedes_albopictus.fna
# Aedes mascarensis (GenBank, chromosome-level)
sed 's/^>/>Amas_|/g' GCA_052575835.1_genomic.fna > Aedes_mascarensis.fna
# Aedes sierrensis (Ochlerotatus)
sed 's/^>/>Asie_|/g' GCA_044231785.1_genomic.fna > Aedes_sierrensis.fna
# Aedes japonicus
sed 's/^>/>Ajap_|/g' GCA_052815935.1_genomic.fna > Aedes_japonicus.fna
# Aedes notoscriptus
sed 's/^>/>Anot_|/g' GCA_040801935.1_genomic.fna > Aedes_notoscriptus.fna
# Aedes koreicus
sed 's/^>/>Akor_|/g' GCA_024533555.2_genomic.fna > Aedes_koreicus.fna


# sed is for simple text substitution. We are only changing the first line information, in this case, the name.
# ^> means "find lines that start with >"

# In the case of genome annotation, since .gtf files are column based, we will use awk in order to modify only the first column that contains the name. 

echo "Concatenating genome files..."

# Clear existing superreference file
> superreference.fna

cat Aedes_aegypti.fna >> superreference.fna
cat Aedes_albopictus.fna >> superreference.fna
cat Aedes_mascarensis.fna >> superreference.fna
cat Aedes_sierrensis.fna >> superreference.fna
cat Aedes_notoscriptus.fna >> superreference.fna
cat Aedes_japonicus.fna >> superreference.fna
cat Aedes_koreicus.fna >> superreference.fna

echo "Superreference files created successfully"

mv superreference.fna /Users/Parsimony/Desktop/jorge_folder/datos/secuencias_Ae_serratus/reference_genomes/aedes_genomes