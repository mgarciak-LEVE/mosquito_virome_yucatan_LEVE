#!/bin/bash

# Script for 
# Author: Jorge Alberto Castro Rodríguez
# Ver. 1.0.0
# 18/05/2026

####==================================####
####          CONFIGURATION           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
${RAW_DATA_DIR}/mapping/${VIRUS_NAME}/${sample_name}/${sample_name}_sorted.bam


####==================================####
####             COVERAGE             ####
####==================================####

