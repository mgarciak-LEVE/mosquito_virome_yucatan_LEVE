# This series of script are designed to download and prepare the databases for taxonomic analysis. 


## Here we will provide the scripts and instructions to download and prepare the databases for viral metagenomics analysis. Command line tools will be used to download, format the databases, and filter the results. 

## Prerequisites.

1. Install the lastest version of DIAMOND (https://github.com/bbuchfink/diamond)
2. Install BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
3. Install wget (https://www.gnu.org/software/wget/)
4. Install pigz (https://zlib.net/pigz/)
5. Install conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Scripts' structure

The scripts are organized into the following main sections:

1. **Database Downloading**: Scripts to download the necessary databases from NCBI and other sources.
2. **Database Preparation**: Scripts to format and prepare the downloaded databases for use with DIAMOND and BLAST.
   - blast_output.sh
   - diamond_output.sh


