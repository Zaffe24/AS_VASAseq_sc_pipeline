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
Firstly, raw FASTQs are demultiplexed into single-cell files, then reads 3'-end are trimmed as per the original [VASA-seq workflow](https://github.com/hemberg-lab/VASAseq_2022/tree/main/I_Gene_expression/a_Mapping). The depletion of ribosomal RNA is performed by mapping reads to a reference FASTA file containing both mouse and human rRNA genes \(`rRNA_ref/rRNA_human_mouse.fa`\). The last step consists in discarding mouse reads that contaminated the human T-ALL sample.

### Installation 
Move to your working directory and clone this repository (change the part between angle brackets):

```shell
git clone https://github.com/Zaffe24/AS_VASAseq_sc_pipeline.git
cd pre_process/ # move to the subfolder for this first module
```
> [!NOTE]
> Run the pipeline from the `pre_process` subfolder and set it as your working directory in `profile/variables.yaml`.

Then create a conda environment named `pre_process` using as template `pre_process.yaml`:

```shell
conda env create -f envs/pre_process.yaml -p </path/to/conda>/pre_process
conda activate pre_process
```
### Execution
Before running the pipeline edit `profile/variables.yaml`. This file allows the user to specify the inputs for the pipeline. Each parameter provides a thorough explanation of its purpose. Moreover, `profile/config.yaml` contains parameters specific to Snakemake and how job submission to the cluster is handled.  

> [!IMPORTANT]
> The ribo-depletion step requires indexing the concatenated human-mouse (hg38+mm38) genome. To build the index, execute `scripts/BBsplit.sh` manually (`BBMap` module required).

Then run:
```shell
snakemake --profile profile/ filtering_cleaned
```

### Output
- `logs/<patient_id>/` contains log reports of each job run on the cluster, organized by rule.
- All the output files are stored in the `<patient_id>/cleaned` directory:
  - List of demultiplexed and mouse-depleted FASTQs.
  - **`units.tsv`** contains the label prefix of all FASTQs generated (`sample` field) and their location (`fq1`). This and the following .tsv files will be needed for modules [III](#iii-micro-exon-discovery) and [IV](#iv-as-event-quantification).
  - **`samples.tsv`** is also required for the following modules, it contains the label prefix of each FASTQ (`sample`) and the respective cluster of origin (`condition`).
  - **`locals.tsv`** contains the path to each FASTQ (`path`), their label prefix (`sample`), and the cluster label to which that cell belongs (`cluster`).
  - **`fractions.tsv`** contains the following order: label prefix for each FASTQ, Proportion of reads mapping to the mouse genome, and proportion of ambiguous reads (mapping to both human and mouse genomes).
  - **`full_metadata.tsv`** merges the information relative to each FASTQ defined in the previous tables.
  - **`Mouse_reads.png`** and **`Ambigous_reads.png`** show boxplots for the proportion of mouse and ambiguous reads in processed and discarded (if a Thr < 1 was set in `variables.yaml`) cells. 
 
  
## II. Transcriptome extension
This module allows the user to extend the (human) genome annotation file (GTF) with novel transcripts inferred directly from the patient/sample data. The single-cell FASTQs are collapsed into pseudo-bulks corresponding to the cell clusters previously defined, and novel transcript annotations are generated for each one of them. The newly annotated transcripts are then appended to the original GTF file that will be used to quantify AS events. For this step, we relied on the pipeline from the original VASA-seq publication.

### Installation
For the installation of this module, please refer to the [original documentation](https://github.com/hemberg-lab/VASAseq_2022/tree/main/II_Alternative_splicing/a_Transcriptome_assembly).

> [!NOTE]
> To avoid confusion between Snakemake configuration files, we suggest installing this pipeline in a different directory.

We ran the pipeline from a dedicated conda environment whose template can be found in [`transcriptome_extension/`](tran scriptome_extension/):

```bash
conda env create -f transcriptome_extension/module_2.yaml -p </path/to/conda>/module_2.yaml

```

## III. Micro-exon discovery

## IV. AS event quantification
