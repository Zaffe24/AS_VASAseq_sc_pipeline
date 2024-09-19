# Alternative Splicing pipeline for VASA-seq
Repository describing the computational workflow employed for the single-cell Alternative Splicing (AS) analysis from [Costea et al., 2024](https://doi.org/10.1101/2024.06.24.600391)

## Overview
Patient-derived xenograft (PDX) samples from pediatric T-cell acute lymphoblastic leukemia (T-ALL) patients were sequenced via VASA-seq, a plate-based total-transcriptome (single-end) scRNA-seq method developed by [Salmen et al., 2022](https://www.nature.com/articles/s41587-022-01361-8). This repository contains the workflow adapted from the original VASA-seq [GitHub repository](https://github.com/hemberg-lab/VASAseq_2022/tree/main/II_Alternative_splicing) for conducting Alternative Splicing (AS) analysis, starting from the raw sequencing output. The whole pipeline is organized into four independent modules:

- [I. Pre-processing](#i-pre-processing)
- [II. Transcriptome extension](#ii-transcriptome-extension)
- [III. Micro-exon discovery](#iii-micro-exon-discovery)
- [IV. AS event quantification](#iv-as-event-quantification)

## I. Pre-processing
The workflow is designed to be run **individually** for each patient/sample via a customized `Snakemake` pipeline.
Firstly, raw FASTQs are demultiplexed into single-cell files, then reads 3'-end are trimmed as per the original [VASA-seq workflow](https://github.com/hemberg-lab/VASAseq_2022/tree/main/I_Gene_expression/a_Mapping). The depletion of ribosomal RNA is performed by mapping reads to a reference FASTA file containing both mouse and human rRNA genes \(`XXXX.fasta`\). The last step consists in discarding mouse reads that contaminated the human T-ALL sample.

### Installation 
Create a conda environment named `pre_process` using as template `pre_process.yaml`:

```shell
conda env create --file=pre_process.yaml --prefix=path/to/conda/pre_process
```

## II. Transcriptome extension

## III. Micro-exon discovery

## IV. AS event quantification
