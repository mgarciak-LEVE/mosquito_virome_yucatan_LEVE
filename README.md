# Mosquito Virome Pipeline

## General Objective.
Characterize viral diversity in Ae. serratus and Ae. taeniorhynchus mosquitoes from the Mexican Republic collected in the Yucatán Peninsula and contrasting environments (conserved, diversified, and urban regions) using total RNA and small RNAs.

### Specific Objectives.

- Characterize the virome of mosquito species collected in the Yucatán Peninsula from 4 habitat types and perform genetic and virological characterization of detected viruses.
- Identify and analyze possible novel viral sequences.
- Determine if there are possible patterns and biogeographic, ecological, and/or taxonomic factors associated with the composition of their viromes.

## Directory Structure

### NFS Storage (Git Repository)

**Location:** `~/git_repos/mosquito_virome_yucatan_LEVE/`

This directory contains **only code, small files, and documentation**. All large data files are stored on Lustre.

```
mosquito_virome_yucatan_LEVE
│
├── containers/ # Container images (Apptainer/Singularity)
│
├── docs/ # Documentation and reports
│ ├── aedes_genomes_specs/ # Aedes genome specifications
│ ├── czid_sequence_reports/ # Final CZID analysis reports
│ ├── czid_sequence_reports_provisional/ # Preliminary reports
│ └── md5_files/ # MD5 checksums for data integrity
│
├── scripts/ # All pipeline scripts
│ ├── aedes_reference_genomes/ # Reference genome indexing
│ ├── czid_analysis/ # CZID data analysis
│ ├── databases/ # Database management scripts
│ ├── individual_analyses/ # Standalone analysis scripts
│ └── pipeline_whole/ # Complete pipeline integration
│
│── project_data / SYMLINK to Lustre directory
│
│── raw_data/ # SYMLINK to Lustre raw data directory
│
└── cat_files/ # SYMLINK to concatenated raw files directory
```

### Lustre

**Location:** `/lustre/scratch126/tol/teams/lawniczak/users/jr46/projects/mosquito_virome_yucatan_LEVE/`

This directory contains **all input data, intermediate files, logs, and results** generated during pipeline execution.

```
mosquito_virome_yucatan_LEVE/
│
├── data/ # All input data
│ ├── raw/
│ │ ├── metadata/ # Sample metadata
│ │ ├── small_RNA/ # Small RNA sequencing data
│ │ └── total_RNA/ # Total RNA sequencing data
│ │   └── cat_files/ # Concatenated files
│ │
│ └── references/ # Reference genomes and databases
│ ├── aedes_super_index/ # Aedes supergenome index
│ ├── databases/
│ │ ├── BLAST/ # BLAST databases
│ │ └── DIAMOND/ # DIAMOND databases
│ └── mosquito_genomes/ # Mosquito genome files
│ └── genomic_files/
│
├── logs/ # Job execution logs
│ ├── assembly/
│ ├── blast/
│ ├── mapping/
│ └── trimming/
│
└── results/ # Analysis results
  ├── aligned/ # Alignment files (BAM/SAM)
  ├── assembly/ # Assembly results
  │ ├── MEGAhit/
  │ ├── metaSPAdes/
  │ ├── rnaSPAdes/
  │ └── statistics/
  ├── blast/ # BLAST results
  ├── czid/ # CZID analysis outputs
  │ └── coverage/ # Coverage statistics
  ├── orf_predicition/ # ORF prediction results
  │ └── czid/
  ├── trimmed/ # Trimmed sequences
  ├── trimmed_qc/ # QC after trimming
  │ ├── fastqc/
  │ ├── multiqc/ 
  │ └── stats/
  └── untrimmed_qc/ # QC before trimming
    ├── fastqc/
    ├── multiqc/ 
    └── stats/
```

