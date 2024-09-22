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
  - **`units.tsv`** contains the label prefix of all FASTQs generated (`sample` field) and their location (`fq1`). This and the following .tsv files will be needed for modules [II](#ii-transcriptome-extension), [III](#iii-micro-exon-discovery) and [IV](#iv-as-event-quantification).
  - **`samples.tsv`** is also required for the following modules, it contains the label prefix of each FASTQ (`sample`) and the respective cluster of origin (`condition`).
  - **`locals.tsv`** contains the path to each FASTQ (`path`), their label prefix (`sample`), and the cluster label to which that cell belongs (`cluster`).
  - **`fractions.tsv`** contains the following order: label prefix for each FASTQ, Proportion of reads mapping to the mouse genome, and proportion of ambiguous reads (mapping to both human and mouse genomes).
  - **`full_metadata.tsv`** merges the information relative to each FASTQ defined in the previous tables.
  - **`Mouse_reads.png`** and **`Ambigous_reads.png`** show boxplots for the proportion of mouse and ambiguous reads in processed and discarded (if a Thr < 1 was set in `variables.yaml`) cells. 
 
  
## II. Transcriptome extension
This module allows the user to extend the (human) genome annotation file (GTF) with novel transcripts inferred directly from the patient/sample data. The single-cell FASTQs are collapsed into pseudo-bulks corresponding to the cell clusters previously defined, and novel transcript annotations are generated for each one of them. The newly annotated transcripts are then appended to the original GTF file that will be used to quantify AS events. For this step, we relied on the pipeline from the original [VASA-seq publication](https://www.nature.com/articles/s41587-022-01361-8).

### Set-up
For the installation of this module, please refer to the [original documentation](https://github.com/hemberg-lab/VASAseq_2022/tree/main/II_Alternative_splicing/a_Transcriptome_assembly).

> [!NOTE]
> To avoid confusion between Snakemake configuration files, we suggest creating a directory specific to this module.

We ran the pipeline from a dedicated conda environment whose template can be found in [`transcriptome_extension/`](transcriptome_extension/):

```bash
conda env create -f transcriptome_extension/module_2.yaml -p </path/to/conda>/module_2.yaml
conda activate module_2
```

### Execution
As before, the parameters in `variables.yaml` and `config.yaml` should be defined according to the user specifics. Our customized templates can be found [here](transcriptome_extension/profile/), otherwise the templates from the original [GitHub repository](https://github.com/hemberg-lab/VASAseq_2022/tree/main/II_Alternative_splicing/a_Transcriptome_assembly/Build) work fine.

> [!IMPORTANT]
> We adapted the original `Snakefile` to satisfy our requirements. Replace the original with [`transcriptome_extension/Snakefile`](transcriptome_extension/Snakefile) in your working directory.

To run the pipeline:

```bash
snakemake --profile profile/ <patient_id>/gffcompare/extended_ref_annotation.gtf
```

### Output
The only relevant output for the next steps is the original GTF file expanded with novel transcript annotations: `<patient_id>/gffcompare/extended_ref_annotation.gtf`.


## III. Micro-exon discovery
This step allows the user to further expand the transcriptome annotation by inferring novel micro-exons (length < 30nt) directly from the data.

### Set-up
We suggest consulting the [original documentation](https://github.com/hemberg-lab/VASAseq_2022/tree/main/II_Alternative_splicing/b_Microexon_annotation) for further information about this step.
We created another dedicated conda environment:

```shell
conda env create -f microexon_discovery/module_3.yaml -p </path/to/conda>/module_3
conda activate module_3
```

### Execution
This and the following modules, rely on the Snakemake pipeline `MicroExonator`. For further information about the inputs required and the rationale behind the workflow, consult the original [MicroExonator page](https://microexonator.readthedocs.io/en/latest/index.html).
We provide our customized configuration files for Snakemake in [microexon_discovery/profile](microexon_discovery/profile).

To run this module with our specifics, set the following parameter in the header of `MicroExonator.smk`:
```python
## header of MicroExonator.smk
configfile : "profile/variables_discovery.yaml"
```
then:
```bash
snakemake -s MicroExonator.smk --profile profile/ Report/out.high_quality.txt
````
### Output
`Report/out.high_quality.txt` lists all the high-quality micro-exons discovered. To append their annotations to the GTF file obtained in the previous module we ran the following script:
```bash
microexon_discovery/add_ME_to_ref.sh "Report/out.high_quality.txt" "</path/module2/patient_id>/gffcompare/extended_ref_annotation.gtf" "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# arg_1 : ouput from module III
# arg 2 : output from module II
# arg 3 : human genome FASTA file
```
This will generate **`Report/out.high_quality_calibrated.gtf`**, the final GTF file incorporating both the annotation of novel transcripts and microexons.

>[!TIP]
>The pipeline will generate all output folders in the current working directory. To avoid the risk of overwriting them when re-running the pipeline on a different patient, we suggest moving the outputs to patient-labeled subfolders after completion.
>
## IV. AS event quantification
In this step, we perform pairwise comparisons across cell pseudo-bulks to compute the list of AS events statistically significant. Moreover, we can compute the percentage-spliced-in (PSI) value for each splicing node at the single-cell level. For further information regarding these outputs, consult the original [pipeline repository](https://github.com/hemberg-lab/VASAseq_2022/tree/main/II_Alternative_splicing/c_AS_quantification).

### Pairwise comparisons
Similarly to the previous module, we defined [`variables_pairwise_comparisons.yaml`](AS_quantification/profile/) as config file for `MicroExonator.smk`:
```python
## MicroExonator.smk header
configfile : "profile/variables_pairwise_comparisons.yaml"
```
Having activated the conda environment from the previous module, run the following code:
```shell
snakemake -s MicroExonator.smk --profile profile/ snakepool
```
#### Output
As output, a .tsv file for each pairwise comparison will be generated in the subfolder **`Whippet/Delta/Single_Cell/Single_nodes/`**. These tables contain information for all AS events detected, regardless if they are statistically significant or not.

### Single-cell PSI
To compute PSI values at the single-cell level instead, define [`variables_psi_matrix.yaml`](AS_quantification/profile/) as config file for  `MicroExonator.smk`:
```python
## MicroExonator.smk header
configfile : "profile/variables_psi_matrix.yaml"
```

Within the same conda environment, run the following code:
```shell
snakemake -s MicroExonator.smk --profile profile/ quant_unpool_single_cell
```

#### Output
This will generate for each cell a .tsv file with PSI values for each splicing node stored in **`Whippet/Quant/Single_Cell/Unpooled`**.

> [!TIP]
> Again, after completion of this module we suggest moving all the output files to a patient/sample-labeled subfolder.
