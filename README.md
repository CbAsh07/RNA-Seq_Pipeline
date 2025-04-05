# RNA-Seq_Pipeline
RNA-Seq pipeline was developed for transcriptomic and genomic analysis

Pipeline Workflow
1. Download Samples &/ store in samples directory

2. Quality Control with FastQC
o Output stored in output/fastqc (Checks base quality, adapter contamination, GC content, etc.)

3. Read Trimming with Trimmomatic
o Adapters removed using TruSeq3-PE-2.fa. (output/trimmed)

4. Alignment with HISAT2
o Reads aligned to the reference genome (GRCh38). (output/alignment)

5. Quantification with featureCounts
o Generates gene counts matrix. (output/featureCounts)

6. Differential Expression Analysis with DESeq2
o Input: featureCounts output. (output/deseq2)

7. Transcript Assembly with StringTie
o Per-sample transcript assembly for gene level analysis. (output/stringtie)

8. Transcript-Level Analysis with Ballgown
o Uses StringTie output to analyze differential expression at transcript level. (output/ballgown)
