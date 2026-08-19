#!/usr/bin/env python3
# Convert ICTV VMR to GRAViTy format and create subset

import pandas as pd
import sys
import glob
import os

# Set directories
OUTPUT_DIR = os.path.expanduser("~/databases/gravity_db")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Find the VMR Excel file in the database directory
vmr_files = glob.glob(os.path.join(OUTPUT_DIR, "*.xlsx")) + glob.glob(os.path.join(OUTPUT_DIR, "*.xls"))

if not vmr_files:
    print("ERROR: No Excel file found in ~/databases/gravity_db/")
    print(f"Please place the VMR Excel file in: {OUTPUT_DIR}")
    sys.exit(1)

vmr_file = vmr_files[0]
print(f"Using file: {vmr_file}")
print(f"Output directory: {OUTPUT_DIR}")

# Read the Excel file
df = pd.read_excel(vmr_file, sheet_name=0)
print(f"Total entries in VMR: {len(df)}")
print("Available columns:", df.columns.tolist())

# Create GRAViTy-compatible table
gravity_df = pd.DataFrame()

# Map columns based on your VMR structure
gravity_df['Order'] = df['Order'].fillna('') if 'Order' in df.columns else ''
gravity_df['Family'] = df['Family'].fillna('') if 'Family' in df.columns else ''
gravity_df['Subfamily'] = df['Subfamily'].fillna('') if 'Subfamily' in df.columns else ''
gravity_df['Genus'] = df['Genus'].fillna('') if 'Genus' in df.columns else ''
gravity_df['Virus name (s)'] = df['Virus name(s)'].fillna('') if 'Virus name(s)' in df.columns else ''
gravity_df['Virus GENBANK accession'] = df['Virus GENBANK accession'].fillna('') if 'Virus GENBANK accession' in df.columns else ''

# Determine Baltimore Group from Genome type
if 'Genome' in df.columns:
    baltimore_map = {
        'dsDNA': 'I',
        'ssDNA': 'II',
        'dsRNA': 'III',
        'ssRNA(+)': 'IV',
        'ssRNA(-)': 'V',
        'ssRNA(RT)': 'VI',
        'dsDNA(RT)': 'VII'
    }
    gravity_df['Baltimore Group'] = df['Genome'].map(baltimore_map).fillna('')
else:
    gravity_df['Baltimore Group'] = ''

# Add other required columns
gravity_df['Virus sequence complete'] = df['Genome coverage'].fillna('Complete') if 'Genome coverage' in df.columns else 'Complete'
gravity_df['Genetic code table'] = '1'
gravity_df['Taxonomic grouping'] = gravity_df['Family'].fillna('Unknown')

# Filter to only include rows with GenBank accessions
gravity_df = gravity_df[gravity_df['Virus GENBANK accession'].notna()]
gravity_df = gravity_df[gravity_df['Virus GENBANK accession'] != '']

# Clean accessions
gravity_df['Virus GENBANK accession'] = gravity_df['Virus GENBANK accession'].str.split('[,;]').str[0].str.strip()

# Save files to output directory
full_output = os.path.join(OUTPUT_DIR, 'virus_description_table.txt')
gravity_df.to_csv(full_output, sep='\t', index=False)
print(f"Created {full_output} with {len(gravity_df)} entries")

# Create subset
print("\nCreating representative subset...")
representatives = []
for family, group in gravity_df.groupby('Family'):
    if family and family != 'Unknown' and family != '':
        complete = group[group['Virus sequence complete'] == 'Complete genome']
        if len(complete) > 0:
            selected = complete.head(2)
        else:
            selected = group.head(2)
        representatives.append(selected)

if representatives:
    subset_df = pd.concat(representatives)
    subset_output = os.path.join(OUTPUT_DIR, 'virus_description_table_subset.txt')
    subset_df.to_csv(subset_output, sep='\t', index=False)
    print(f"Created {subset_output} with {len(subset_df)} viruses")
    print(f"Representing {len(subset_df['Family'].unique())} families")
else:
    print("WARNING: No representatives found, using full table")
    subset_output = os.path.join(OUTPUT_DIR, 'virus_description_table_subset.txt')
    gravity_df.to_csv(subset_output, sep='\t', index=False)
