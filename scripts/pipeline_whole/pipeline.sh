#!/bin/bash

SECONDS=0

# If run on a different computer, make sure the Working Directory path is correct and the permissions to run the code are correct as well
# First of all, we establish a script directory

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" # We assume that the scripts are located in a subdirectory "scripts/"
# Works regardless of the current working directory.
PROJECT_ROOT="$SCRIPT_DIR"

# RELATIVE PATHS
RAW_DATA_DIR="${PROJECT_ROOT}/data/Aedes_RNA_sequences"

# RESULTS DIRECTORIES
RESULTS_DIR="${PROJECT_ROOT}/results/aedes_analysis"
UNTRIMMED_QC_DIR="${RESULTS_DIR}/untrimmed"
TRIMMED_DIR="${RESULTS_DIR}/trimmed"
TRIMMED_QC_DIR="${RESULTS_DIR}/trimmed_qc"
ALIGNED_DIR="${RESULTS_DIR}/aligned"
ASSEMBLY_DIR="${RESULTS_DIR}/assembly"
BLAST_DIR="${RESULTS_DIR}/blast"

# DATABASES AND REFERENCE GENOMES
REFERENCE_DIR="${PROJECT_ROOT}/reference_genomes/aedes_genomes/aedes_super_index"
BLAST_DB="${PROJECT_ROOT}/databases/blast/ref_viruses_rep_genomes"
DIAMOND_DB="${PROJECT_ROOT}/databases/diamond/viral_proteins.dmnd"

# LOGS DIRECTORY
LOGS_DIR="${PROJECT_ROOT}/logs"

# WORKING AND OUTPUT DIRECTORIES
WORKDIR="$RAW_DATA_DIR"
OUTDIR="$RESULTS_DIR"


threads=16 # Thread definition

# A pipline was used to process some Aedes serratus sequences

# SAMPLE NAME EXTRACTION
sample_names() {
    local filename=$1                    # Takes the local file name with the current format "PM1234_SX_R1" and its file path
    local option="${2:-sample}"          # We just need the sample name
    local basename=$(basename "$filename" .fastq)  # Removes basename (directory path) and .fastq 
    basename=$(echo "$basename" | sed 's/_[0-9]\{3\}$//')

    # Name extraction 
    local sample_part=$(echo "$basename" | cut -d'_' -f1)  # Gets "PMxxxx"
    local read_part=$(echo "$basename" | grep -o 'R[12]$')  # Gets "R1" or "R2"
    
    # Gives back 4 different parts
    case $option in
        "sample") echo "$sample_part" ;;     # Just "PM2469"
        "read") echo "$read_part" ;;         # Just "R1" or "R2"  
        "full") echo "${sample_part}_${read_part}" ;;  # "PM2469_R1"
        "with_lane") echo "$basename" ;;     # "PM2469_S1_R1"
    esac
}

# FASTQC FUNCTION
run_fastqc() {
    # Runs FastQC based on the files that have a .fastq format
    local input_dir=$1 # Takes the first argument
    local output_dir=$2 # Takes the second argument

    echo "Running FastQC on files in $input_dir..."
    mkdir -p "${output_dir}"

    find "$input_dir" -name "*.fastq" | while read -r file; do # Searches the input directory for the files ending with .fastq and meanwhile... \\ 
        echo "Processing $file..."
        fastqc "$file" -o "$output_dir" # Runs fastqc on the current file ($file) and saves the results on output_dir
    done
}


# MULTIQC FUNCTION
run_multiqc() {
    #Runs MultiQC
    local input_dir=$1 # Takes the first argument
    local output_dir=$2 # Takes the second argument

    conda init
    echo "Running MultiQC on $input_dir..."
    conda activate multiqc_env # MultiQC conda environment 
    multiqc "$input_dir" -o "$output_dir"
    conda deactivate
}

# TRIMMING FUNCTION
run_trimming() {
    #Runs trimming for paired-end reads. It distinguishes between the R1 and R2 reads
    echo "Starting trimming process..."
    mkdir -p "${TRIMMED_DIR}" # -p flag means to create parent directories if they dont exist 

    for R1_file in "$WORKDIR"/*_R1.fastq; do   # Does a loop. Finds every file in the working directory that matches the pattern "*_R1.fastq" and assigns each file path to the variable R1
        R2_file="${R1_file/_R1.fastq/_R2.fastq}" # Parameter expansion: ${variable/pattern/replacement}
    
    newR1=$(sample_names "$R1_file" "full")
    newR2=$(sample_names "$R2_file" "full")


    conda activate trimmomatic 
    echo "Input files: $(basename "$R1_file") and $(basename "$R2_file")"
    echo "Processing: $newR1 and $newR2"
    adapters="/Users/Parsimony/miniconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa" 
    
    # For the forward and reverse files ($R1 and $R2)
    trimmomatic PE -threads 10 -phred33 \
        "$R1_file" "$R2_file" \
        "${TRIMMED_DIR}/${newR1}_paired.fastq" \
        "${TRIMMED_DIR}/${newR1}_unpaired.fastq" \
        "${TRIMMED_DIR}/${newR2}_paired.fastq" \
        "${TRIMMED_DIR}/${newR2}_unpaired.fastq" \
        ILLUMINACLIP:${adapters}:2:30:10:2:keepBothReads \
        LEADING:10 TRAILING:10 SLIDINGWINDOW:3:25 MINLEN:50 # ILUMMINACLIP is for adapter clipping
    done   
    conda deactivate
}

# MAPPING FUNCTION
run_mapping() {
    echo "Starting alignment process..."
    mkdir -p "${ALIGNED_DIR}" # Output directory for the aligned sequences
    REFERENCE="$REFERENCE_DIR"

    ulimit -n 4096

    for R1 in "$OUTDIR"/trimmed/*_R1_paired.fastq; do
        # Get base name
        R2="${R1/_R1_paired.fastq/_R2_paired.fastq}"

    sample=$(sample_names "$R1" "sample")


    echo "Using STAR index from: $REFERENCE"
    echo "Aligning $sample..." # Alignment done by STAR 
    echo "Files: $(basename "$R1") and $(basename "$R2")"

    conda activate star_env # STAR conda environment

    STAR \
        --runMode alignReads \
        --genomeDir "$REFERENCE" \
        --readFilesIn "$R1" "$R2" \
        --outFileNamePrefix "${ALIGNED_DIR}/${sample}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outReadsUnmapped Fastx \
        --quantMode GeneCounts \
        --outFilterMismatchNoverLmax 0.1 \
        --runThreadN $threads \
        --limitBAMsortRAM 30000000000 \
        
    conda deactivate
    done
}

# ASSEMBLY INFORMATION
mapping_stats() {
    echo "Generating mapping statistics..."
    mkdir -p "${ALIGNED_DIR}/statistics"
    
    #Creates CSV file
    stats_file="${OUTDIR}/aligned/statistics/mapping_summary.csv"
    #Specifies the file column headers
    echo "Sample,Input_Reads,Average_Length,Uniquely_Mapped_Reads,Uniquely_Percent,Total_Multimapped,Total_Multimapped_Percent,Unmapped_Reads,Unmapped_Percent,Unmapped_Too_Short,Unmapped_Other" > "$stats_file"


    for star_log in "$OUTDIR/aligned"/*Log.final.out; do
        if [[ -f "$star_log" ]]; then
            sample=$(basename "$star_log" | sed 's/_Log.final.out//')

                    # Parse statistics
                    input_reads=$(grep "Number of input reads" $star_log | awk '{print $6}')
                    avg_length=$(grep "Average input read length" $star_log | awk '{print $6}')
                    #Uniquely mapped
                    uniquely_mapped_reads=$(grep "Uniquely mapped reads number" $star_log  | awk '{print $6}')
                    uniquely_percent=$(grep "Uniquely mapped reads %" $star_log  | awk '{print $6}' | sed 's/%//')
                    #Multimapped
                    multi_mapped1=$(grep "Number of reads mapped to multiple loci" $star_log  | awk '{print $9}')
                    multi_mapped1_percent=$(grep "% of reads mapped to multiple loci" $star_log  | awk '{print $9}' | sed 's/%//')
                    multi_mapped2=$(grep "Number of reads mapped to too many loci" $star_log  | awk '{print $10}')
                    multi_mapped2_percent=$(grep "% of reads mapped to too many loci" $star_log  | awk '{print $10}' | sed 's/%//')
                    total_multi_mapped=$(($multi_mapped1 + $multi_mapped2))
                    total_multi_percent=$(echo "scale=2; $total_multi_mapped * 100 / $input_reads" | bc)

                    #Unmapped
                    unmapped_mismatches=$(grep "Number of reads unmapped: too many mismatches" $star_log  | awk '{print $9}')
                    unmapped_short=$(grep "Number of reads unmapped: too short" $star_log  | awk '{print $8}')
                    unmapped_other=$(grep "Number of reads unmapped: other" $star_log  | awk '{print $7}')
                    total_unmapped=$(($unmapped_mismatches + $unmapped_short + $unmapped_other))
                    unmapped_percent=$(echo "scale=2; $total_unmapped * 100 / $input_reads" | bc)

                    echo "$sample,$input_reads,$avg_length,$uniquely_mapped_reads,$uniquely_percent,$total_multi_mapped,$total_multi_percent,$total_unmapped,$unmapped_percent,$unmapped_short,$unmapped_other" >> "$stats_file"
                fi
    done
    
    echo "Statistics saved to: $stats_file"
    echo "=== Mapping Statistics ==="
    column -t -s, "$stats_file" 

}



# ASSEMBLY FUNCTION 
run_assembly() {
    echo "Starting assembly process..."
    mkdir -p "${OUTDIR}/assembly" # Output directory for assembly
    mkdir -p "${OUTDIR}/assembly/fastq"

    for R1 in "$OUTDIR"/aligned/*_Unmapped.out.mate1; do
        R2="${R1/_Unmapped.out.mate1/_Unmapped.out.mate2}"
        
        sample=$(basename $R1 | cut -d'_' -f1)

        # Sample specific output directories
        mkdir -p "$OUTDIR/assembly/rnaSPAdes/$sample"
        mkdir -p "$OUTDIR/assembly/metaSPAdes/$sample"
        echo "Cleaning up previous MEGAhit results for $sample..."
        rm -rf "$OUTDIR/assembly/MEGAhit/$sample" 2>/dev/null
        

        R1_fastq="$OUTDIR/assembly/fastq/${sample}_R1.fastq"
        R2_fastq="$OUTDIR/assembly/fastq/${sample}_R2.fastq"

        echo "Extracting unmapped reads from BAM format to FASTQ"
        
        conda activate samtools_env
        samtools fastq "$R1" > "$R1_fastq"
        samtools fastq "$R2" > "$R2_fastq"
        conda deactivate

        # rnaSPAdes --------
        # rnaSPAdes is good for transcriptome analyses.
        echo "Starting assembly through rnaSPAdes..."
        # Input data as every corresponding unmapped read file.

        echo "Starting rnaSPAdes assembly for sample $sample..."
        conda activate SPADES_env
        rnaspades.py \
        -1 "$R1_fastq" \
        -2 "$R2_fastq" \
        -o "$OUTDIR"/assembly/rnaSPAdes/$sample \
        -t $threads
        conda deactivate
        echo "rnaSPAdes assembly finished for sample $sample"

        # metaSPAdes. --------
        # metaSPAdes is good for small genomes.
        echo "Starting assembly through metaSPAdes..." 

        # Input files were the same unmapped reads with the same format. 
        echo "Starting metaSPAdes assembly for sample $sample..."
        conda activate SPADES_env
        metaspades.py \
        -1 "$R1_fastq" \
        -2 "$R2_fastq" \
        -o "$OUTDIR"/assembly/metaSPAdes/$sample \
        -t $threads
        conda deactivate
        echo "metaSPAdes assembly finished for sample $sample"

        # MEGAhit. --------
        echo "Starting assembly through MEGAhit..."

        echo "Starting MEGAhit assembly for sample $sample..."
        conda activate MEGAhit_env 
        megahit \
        -1 "$R1_fastq" \
        -2 "$R2_fastq" \
        -o "$OUTDIR"/assembly/MEGAhit/$sample \
        -t 4 \
        -m 30
        conda deactivate
        echo "MEGAhit assembly finished for sample $sample"

    done

    echo "All assembly processes completed"
}

# ASSEMBLY INFORMATION
assembly_stats() {
    echo "Generating assembly statistics with seqkit..."
    mkdir -p "${OUTDIR}/assembly/statistics"
    
    #Creates CSV file
    stats_file="${OUTDIR}/assembly/statistics/assembly_summary_seqkit.csv"
    #Specifies the file column headers
    echo "Sample,Assembler,File,Num_Contigs,Total_Length,Max_Length,Min_Length,Avg_Length,N50,GC_Percent" > "$stats_file"

    conda activate seqkit_env

    # rnaSPAdes statistics
    for sample_dir in "$OUTDIR/assembly/rnaSPAdes"/*/; do
        if [[ -d "$sample_dir" ]]; then # If its a directory...
            sample=$(basename "$sample_dir")
            
            # RNAspAdes typically produces these files:
            for contig_file in "$sample_dir/transcripts.fasta" \
                              "$sample_dir/soft_filtered_transcripts.fasta" \
                              "$sample_dir/hard_filtered_transcripts.fasta"; do
                if [[ -f "$contig_file" ]]; then # If the file exists...
                    echo "Processing $(basename "$contig_file") for $sample - rnaSPAdes"
                    
                    # Get statistics with seqkit
                    stats=$(seqkit stats "$contig_file" -T | awk 'NR==2')  # -T seqkit data output in tabular format. NR==2 means that it only gets the data row
                    
                    # Parse statistics
                    num_contigs=$(echo "$stats" | awk '{print $4}')
                    total_length=$(echo "$stats" | awk '{print $5}')
                    min_length=$(echo "$stats" | awk '{print $6}')
                    max_length=$(echo "$stats" | awk '{print $7}')
                    avg_length=$(echo "$stats" | awk '{print $8}')

                    # N50 calculation
                    # seqkit fx2tab -n -l converts FASTA to table with names and lengths.
                    n50=$(seqkit fx2tab -n -l "$contig_file" | awk '
                    {len[NR]=$2; sum+=$2} 
                    END {
                        half=sum/2;
                        for(i=1;i<=NR;i++) {
                            running_sum+=len[i];
                            if(running_sum>=half) {
                                print len[i];
                                exit;
                            }
                        }
                    }')
                    #len[NR]=$2; sum+=$2: Stores lengths in array, calculates total sum
                    #half=sum/2: Calculates halfway point
                    #running_sum+=len[i]: Adds lengths sequentially
                    #if(running_sum>=half): Checks when we cross the halfway point
                    #print len[i]: Outputs the N50 value


                    # GC content
                    #seqkit fx2tab -n -g gets GC percentage for each sequence
                    gc_percent=$(seqkit fx2tab -n -g "$contig_file" | awk '
                    {sum_gc+=$2; count++} 
                    END {if(count>0) printf "%.2f", sum_gc/count; else print "0"}')

                    filename=$(basename "$contig_file")
                    echo "$sample,rnaSPAdes,$filename,$num_contigs,$total_length,$max_length,$min_length,$avg_length,$n50,$gc_percent" >> "$stats_file"
                fi
            done
        fi
    done


    # metaSPAdes statistics

    for sample_dir in "$OUTDIR/assembly/metaSPAdes"/*/; do
        if [[ -d "$sample_dir" ]]; then
            sample=$(basename "$sample_dir")
            
            # metaSPAdes typically has these files
            for contig_file in "$sample_dir/contigs.fasta" "$sample_dir/scaffolds.fasta"; do
                if [[ -f "$contig_file" ]]; then
                    file_desc=$(basename "$contig_file" .fasta)
                    echo "Processing $file_desc for $sample - metaSPAdes"
                    
                    # Get statistics with seqkit
                    stats=$(seqkit stats "$contig_file" -T | awk 'NR==2')
                    
                    # Parse statistics
                    num_contigs=$(echo "$stats" | awk '{print $4}')
                    total_length=$(echo "$stats" | awk '{print $5}')
                    min_length=$(echo "$stats" | awk '{print $6}')
                    max_length=$(echo "$stats" | awk '{print $7}')
                    avg_length=$(echo "$stats" | awk '{print $8}')


                    # N50 calculation
                    n50=$(seqkit fx2tab -n -l "$contig_file" | awk '
                    {len[NR]=$2; sum+=$2} 
                    END {
                        half=sum/2;
                        for(i=1;i<=NR;i++) {
                            running_sum+=len[i];
                            if(running_sum>=half) {
                                print len[i];
                                exit;
                            }
                        }
                    }')
                    
                    # GC content
                    gc_percent=$(seqkit fx2tab -n -g "$contig_file" | awk '
                    {sum_gc+=$2; count++} 
                    END {if(count>0) printf "%.2f", sum_gc/count; else print "0"}')

                    
                    echo "$sample,metaSPAdes,$file_desc,$num_contigs,$total_length,$max_length,$min_length,$avg_length,$n50,$gc_percent" >> "$stats_file"
                fi
            done
        fi
    done
    
    echo "Statistics saved to: $stats_file"
    echo "=== Assembly Statistics ==="
    column -t -s, "$stats_file"

    conda deactivate   

}


ntblast() {
    mkdir -p "${BLAST_DIR}/ntblast"

     echo "Running BLAST on assembly results..."

    conda activate blast_env

    # rnaSPAdes
    blastn -query "${ASSEMBLY_DIR}/rnaSPAdes/hard_filtered_transcripts.fasta" \
        -db "$BLAST_DB" \
        -out "${BLAST_DIR}/ntblast/rnaSPAdes_blast_results.txt" \
        -outfmt 6 -evalue 1e-5 -num_threads 8

    # metaSPAdes
    blastn -query "${ASSEMBLY_DIR}/metaSPAdes/contigs.fasta" \
        -db "$BLAST_DB" \
        -out "${BLAST_DIR}/ntblast/metaSPAdes_blast_results.txt" \
        -outfmt 6 -evalue 1e-5 -num_threads 8

    # MEGAhit
    blastn -query "${ASSEMBLY_DIR}/MEGAhit/final.contigs.fa" \
        -db "$BLAST_DB" \
        -out "${BLAST_DIR}/ntblast/MEGAhit_blast_results.txt" \
        -outfmt 6 -evalue 1e-5 -num_threads 8

    conda deactivate
}

nrblast() {
    mkdir -p "${BLAST_DIR}"
    mkdir -p "${BLAST_DIR}/nrblast"
    

    echo "Running DIAMOND on assembly results..."

    conda activate diamond_env

    # rnaSPAdes
    diamond blastx -q "${OUTDIR}/assembly/rnaSPAdes/hard_filtered_transcripts.fasta" \
            -d "$DIAMOND_DB" \
            --o "${BLAST_DIR}/nrblast/rnaSPAdes_diamond_results.txt" \
            --ultra-sensitive --evalue 1e-5 --threads 8 --outfmt 6

    # metaSPAdes
    diamond blastx -q "${OUTDIR}/assembly/metaSPAdes/contigs.fasta" \
            -d "$DIAMOND_DB" \
            --o "${BLAST_DIR}/nrblast/metaSPAdes_diamond_results.txt" \
            --ultra-sensitive --evalue 1e-5 --threads 8 --outfmt 6

    # MEGAhit
    diamond blastx -q "${OUTDIR}/assembly/MEGAhit/final.contigs.fa" \
            -d "$DIAMOND_DB" \
            --o "${BLAST_DIR}/nrblast/MEGAhit_diamond_results.txt" \
            --ultra-sensitive --evalue 1e-5 --threads 8 --outfmt 6

    conda deactivate
}


# SEQUENCE ANALYSIS PIPELINE <-----------
main() {
    eval "$(/Users/Parsimony/miniconda3/bin/conda shell.bash hook)"
    
# Raw QC analysis
    run_fastqc "$WORKDIR" "$UNTRIMMED_QC_DIR"
    run_multiqc "$UNTRIMMED_QC_DIR" "$RESULTS_DIR"
    
    # Trimming process
    run_trimming
    
    # Post-trimming QC analysis
    run_fastqc "$TRIMMED_DIR" "$TRIMMED_QC_DIR"
    run_multiqc "$TRIMMED_QC_DIR" "$RESULTS_DIR"

    # Mapping process
    run_mapping

    # Mapping statistics
    mapping_stats

    # Assembly process
    run_assembly

    # Assembly statistics
    assembly_stats

    # Viral identification
    ntblast
    nrblast


    # Report time
    duration=$SECONDS
    echo "Pipeline finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
}

main 
