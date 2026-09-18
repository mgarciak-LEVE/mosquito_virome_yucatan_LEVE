#!/usr/bin/env python2.7

import pandas as pd

# Read the Excel file, forcing all columns to string to avoid Excel corruption
df = pd.read_excel("VMR_MSL41.v1.20260729.xlsx", dtype=str)

# Keep ONLY the 9 required columns
required_cols = [
    "Virus GENBANK accession",
    "Realm", "Kingdom", "Phylum", "Class",
    "Order", "Family", "Genus", "Species"
]

df[required_cols].to_csv("VMR_MSL41_vitap.csv", index=False)
