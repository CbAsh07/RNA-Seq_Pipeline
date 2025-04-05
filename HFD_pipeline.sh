# Define paths for input/output directories
SAMPLES_DIR="samples"
mkdir -p "$SAMPLES_DIR"
OUTPUT_DIR="output"
mkdir -p "$OUTPUT_DIR"
# Specify directory where genome index and GTF files are downloaded
GENOME_DIR="hisat2-2.2.1"
# Specify directory for R scripts and sample info
DESEQ2_SCRIPT="DESEq2.R"
BALLGOWN_SCRIPT="run_ballgown.R"
SAMPLE_INFO="sample_info.txt"
#Specify output files
FASTQC_OUTPUT="$OUTPUT_DIR/fastqc"
TRIMMED_DIR="$OUTPUT_DIR/trimmed"
ALIGNMENT_DIR="$OUTPUT_DIR/alignment"
FEATURECOUNTS_OUTPUT="$OUTPUT_DIR/featureCounts"
DESEQ2_OUTPUT="$OUTPUT_DIR/deseq2"
STRINGTIE_OUTPUT="$OUTPUT_DIR/stringtie"
BALLGOWN_OUTPUT="$OUTPUT_DIR/ballgown"
LOG="$OUTPUT_DIR/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

mkdir -p "$FASTQC_OUTPUT" "$TRIMMED_DIR" "$ALIGNMENT_DIR" "$FEATURECOUNTS_OUTPUT" "$DESEQ2_OUTPUT" "$STRINGTIE_OUTPUT" "$BALLGOWN_OUTPUT"

# Download RNA-Seq samples into the samples directory
echo "Downloading RNA-Seq samples from sample_info.txt..."

if [[ ! -f "$SAMPLE_INFO" ]]; then
    echo "Error: $SAMPLE_INFO not found!"
    exit 1
fi

declare -A SAMPLES_SET  # Use associative array to collect unique sample IDs

# Read the sample_info.txt (tab-separated: sample, condition, read1_url, read2_url)
while IFS=$'\t' read -r sample condition read1_url read2_url; do
    if [[ "$sample" == "sample" ]]; then
        continue  # Skip header
    fi

    # Define file paths
    read1_file="$SAMPLES_DIR/$(basename "$read1_url")"
    read2_file="$SAMPLES_DIR/$(basename "$read2_url")"

    if [[ -f "$read1_file" && -f "$read2_file" ]]; then
        echo "Skipping download for $sample (files already exist)"
    else
        echo "Downloading files for $sample..."
        wget -P "$SAMPLES_DIR" "$read1_url"
        wget -P "$SAMPLES_DIR" "$read2_url"
    fi

    # Collect sample names into the associative array
    SAMPLES_SET["$sample"]=1

done < "$SAMPLE_INFO"

# Convert associative array keys into a normal array
SAMPLES=("${!SAMPLES_SET[@]}")

# Show collected samples
echo "Samples detected: ${SAMPLES[*]}"

# Quality Control with FastQC
echo "Running FastQC for quality control..."
for SAMPLE in ${SAMPLES[@]}; do
  if [ ! -f "$FASTQC_OUTPUT/${SAMPLE}_1_fastqc.zip" ] || [ ! -f "$FASTQC_OUTPUT/${SAMPLE}_2_fastqc.zip" ]; then
    fastqc -o "$FASTQC_OUTPUT" "$SAMPLES_DIR/${SAMPLE}_1.fastq.gz" "$SAMPLES_DIR/${SAMPLE}_2.fastq.gz"
  else
    echo "FastQC output already exists for $SAMPLE, skipping..."
  fi
done

# Run Trimmomatic
echo "Running Trimmomatic for trimming..."
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Processing $SAMPLE..."
  
  # Check if output files already exist
  if [ ! -f "$TRIMMED_DIR/${SAMPLE}_1_trimmed.fastq.gz" ] || [ ! -f "$TRIMMED_DIR/${SAMPLE}_2_trimmed.fastq.gz" ]; then
    # Run Trimmomatic
    ADAPTERS="/opt/homebrew/opt/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
    trimmomatic PE -phred33 \
    "$SAMPLES_DIR/${SAMPLE}_1.fastq.gz" "$SAMPLES_DIR/${SAMPLE}_2.fastq.gz" \
    "$TRIMMED_DIR/${SAMPLE}_1_trimmed.fastq.gz" "$TRIMMED_DIR/${SAMPLE}_1_unpaired.fastq.gz" \
    "$TRIMMED_DIR/${SAMPLE}_2_trimmed.fastq.gz" "$TRIMMED_DIR/${SAMPLE}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:27 TRAILING:27 MINLEN:10


    echo "$SAMPLE trimming complete!"
  else
    echo "Trimmomatic output already exists for $SAMPLE, skipping..."
  fi
done

echo "All samples processed!"

# Alignment with HISAT2
echo "Running HISAT2 for alignment..."

for SAMPLE in ${SAMPLES[@]}; do
  if [ ! -f "$ALIGNMENT_DIR/${SAMPLE}_sorted.bam" ]; then
    hisat2 -x "$GENOME_DIR/genome" \
         -1 "$TRIMMED_DIR/${SAMPLE}_1_trimmed.fastq.gz" \
         -2 "$TRIMMED_DIR/${SAMPLE}_2_trimmed.fastq.gz" \
         -S "$ALIGNMENT_DIR/${SAMPLE}_aligned.sam"

    # Convert SAM to BAM and sort with SAMtools
    samtools view -bS "$ALIGNMENT_DIR/${SAMPLE}_aligned.sam" | samtools sort -o "$ALIGNMENT_DIR/${SAMPLE}_sorted.bam"
    samtools index "$ALIGNMENT_DIR/${SAMPLE}_sorted.bam"

    # Optionally remove the SAM file to save space
    rm "$ALIGNMENT_DIR/${SAMPLE}_aligned.sam"
  else
    echo "Alignment already done for $SAMPLE, skipping..."
  fi
done

# Create a variable to store sorted BAM file paths
SORTED_BAM_FILES=""

for SAMPLE in ${SAMPLES[@]}; do
  SORTED_BAM="$ALIGNMENT_DIR/${SAMPLE}_sorted.bam"
  SORTED_BAM_FILES+=" $SORTED_BAM"  # Append each BAM file path to the variable
done

# Quantification with featureCounts
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
COUNTS_FILE="counts_${TIMESTAMP}.txt"

echo "Running featureCounts for quantification..."
if [ ! -f "$FEATURECOUNTS_OUTPUT/$COUNTS_FILE" ]; then
    featureCounts -T 4 -a "$GENOME_DIR/Homo_sapiens.GRCh38.113.chr.gtf.gz" \
                  -o "$FEATURECOUNTS_OUTPUT/$COUNTS_FILE" \
                  $SORTED_BAM_FILES 

    # Clean up column names in counts file
    awk -F'\t' 'NR==2 {
        for (i=7; i<=NF; i++) {
            sub(/^output\/alignment\//, "", $i)
            sub(/_sorted\.bam$/, "", $i)
        }
    } 1' OFS='\t' "$FEATURECOUNTS_OUTPUT/$COUNTS_FILE" > "$FEATURECOUNTS_OUTPUT/counts_fixed_${TIMESTAMP}.txt" \
    && mv "$FEATURECOUNTS_OUTPUT/counts_fixed_${TIMESTAMP}.txt" "$FEATURECOUNTS_OUTPUT/$COUNTS_FILE"
else
    echo "featureCounts output already exists, skipping..."
fi

# Ensure FeatureCounts output exists
if [ ! -f "$FEATURECOUNTS_OUTPUT/$COUNTS_FILE" ]; then
    echo "Error: FeatureCounts output not found at $FEATURECOUNTS_OUTPUT/$COUNTS_FILE"
    exit 1
fi

echo "FeatureCounts completed successfully: $COUNTS_FILE"

# Run DESeq2 R script
echo "Running DESeq2 for differential expression analysis..."
Rscript "$DESEQ2_SCRIPT" "$FEATURECOUNTS_OUTPUT/$COUNTS_FILE" "$DESEQ2_OUTPUT" "$SAMPLE_INFO"

# Check if the analysis completed successfully
if [ $? -eq 0 ]; then
    echo "DESeq2 analysis completed successfully."
    echo "Results can be found in: $DESEQ2_OUTPUT"
else
    echo "Error: DESeq2 analysis failed."
    exit 1
fi

# RUN StringTie (for Ballgown)
echo "Running StringTie for transcript assembly (Ballgown preparation)..."

for SAMPLE in ${SAMPLES[@]}; do
    SAMPLE_DIR="$STRINGTIE_OUTPUT/$SAMPLE"
    OUTPUT_GTF="$SAMPLE_DIR/${SAMPLE}.gtf"

    # Check if the output GTF file already exists
    if [[ -f "$OUTPUT_GTF" ]]; then
        echo "Skipping $SAMPLE: Output already exists at $OUTPUT_GTF"
        continue
    fi

    mkdir -p "$SAMPLE_DIR"

    stringtie "$ALIGNMENT_DIR/${SAMPLE}_sorted.bam" \
        -o "$OUTPUT_GTF" \
        -p 8 -G "$GENOME_DIR/Homo_sapiens.GRCh38.113.chr.gtf" \
        -B  # This creates Ballgown-ready files (ctab files inside each sample folder)

    if [[ $? -eq 0 ]]; then
        echo "StringTie completed for $SAMPLE."
    else
        echo "Error: StringTie failed for $SAMPLE."
        exit 1
    fi
done

# Run Ballgown analysis using the R script
echo "Running Ballgown analysis using run_ballgown.R..."

if [[ "$OSTYPE" == "darwin"* ]]; then
    RSCRIPT="/Library/Frameworks/R.framework/Resources/bin/Rscript"
else
    RSCRIPT="Rscript"
fi

$RSCRIPT "$BALLGOWN_SCRIPT" "$OUTPUT_DIR/ballgown"

if [ $? -eq 0 ]; then
    echo "Ballgown analysis completed successfully."
else
    echo "Error: Ballgown analysis failed."
    exit 1
fi

echo "Both featureCounts-DESeq2 and StringTie-Ballgown workflows completed!"

echo "RNA-Seq pipeline completed."