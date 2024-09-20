  ###########################################################################
 #### Snakemake pipeline for the pre-processing of raw VASA-seq output ####
#########################################################################
# by Pietro Zafferani


#### import all packages ####
#### make sure that the pipeline is using the Python from the conda env ####

import pandas as pd
import random
from collections import defaultdict
import csv
import glob2
import sys

#### import user-defined variables ####
configfile : "profile/variables.yaml"

#### define plate and cell IDs ####
plates=config["plates"].split(';')
plates=[p.strip() for p in plates]
cells=["{0:03}".format(i) for i in range(1,385)]


tempo_out=os.path.join(config["tmp"],config["patient"])

  ###############################################################################
 ############################# DEMULTIPLEXING ##################################
###############################################################################

#### Running this rule on an HPC should take around 1h per plate ####
for plate in plates:
  rule:
    name:plate + '_demultiplexing'
    
    params:
      patient=config["patient"],
      p1 = tempo_out,
      p2 = config["path_to_plates"],
      p3 = config["working_dir"],
      pl = plate,
    
    resources:
      mem_mb=1000
    
    ### demultiplexed FASTQs will be stored in temporary folder ###  
    output:
      expand(os.path.join(tempo_out, "raw_fastq" , plate+"_{cell_id}_cbc.fastq.gz"), cell_id=cells),

    shell:
      """
      echo -e plate: {params.pl} 
      scripts/demultiplexing.sh {params.p1} '{params.p2}' {params.p3} {params.pl}
      """


  ###############################################################################
 ################################## TRIMMING ###################################
###############################################################################

rule trimming:
  input:
    os.path.join(tempo_out,"raw_fastq","{plate}_{cell_id}_cbc.fastq.gz")
      
  params:
    patient=config["patient"],
    trimal=config["trimal"],
    cutadapt=config["cutadapt"],
    outdir=os.path.join(tempo_out,"trimmed"),
      
  resources:
    load=int(config["fast_jobs"])

  output:
    os.path.join(tempo_out,"trimmed/{plate}_{cell_id}_cbc_trimmed_homoATCG.fq.gz"),
    
  shell:
    """
    module load Trim_Galore/0.6.7-GCCcore-10.3.0      
    module load cutadapt/3.4-GCCcore-10.3.0
    
    scripts/trimming.sh {input} {params.outdir} {params.trimal} {params.cutadapt}
    """
    
    
  ###############################################################################
 ########################### RIBOSOME DEPLETION ################################
###############################################################################


rule ribo_depletion:
  input:
    os.path.join(tempo_out, "trimmed/{plate}_{cell_id}_cbc_trimmed_homoATCG.fq.gz")
  
  params:
    patient=config["patient"],
    hs_mm_ribo=config["hs_mm_ribo"],
    bwa=config["bwa"],
    samtools=config["samtools"],
    outdir=os.path.join(tempo_out, "ribo_depleted"),
    
  resources:
    load=int(config["fast_jobs"])
    
  output:
    os.path.join(tempo_out, "ribo_depleted/{plate}_{cell_id}.nonRibo.fastq.gz")
  
  shell:
    """
    pref_out=$(basename {input}) #pass input prefix to output prefix
    echo $pref_out
    module load BWA/0.7.17-GCCcore-11.2.0 
    module load SAMtools/1.13-GCC-11.2.0
    scripts/ribo-bwamem.sh {params.hs_mm_ribo} {input} {params.outdir}/"$pref_out" {params.bwa} {params.samtools} "y" "scripts"
    sleep 2
    
    file_size=$(stat -c%s {output})
    min=100
    if [ $file_size -gt $min ]; then
      echo -e 'file size in bytes:' "$file_size"
      echo ribo-depletion for {input} completed
    else
      echo -e 'file size in bytes:' "$file_size" 'insufficient'
      exit -1
    fi
    """


  ###############################################################################
 ########################### MOUSE READ FILTERING ##############################
###############################################################################

rule read_cleaning:
  input:
    os.path.join(tempo_out, "ribo_depleted/{plate}_{cell_id}.nonRibo.fastq.gz")
    
  params:
    patient=config["patient"],
    outdir=os.path.join(tempo_out,"cleaned"),
    ref_path=config["ref_path"],
    bbmap=config["bbmap"],
    
  resources:
    load=int(config["slow_jobs"])
    
  output:
    os.path.join(tempo_out, 'cleaned/{plate}_{cell_id}.cleaned.fq.gz'),

  shell:
    """
    module load BBMap    ##try without
    scripts/pdx.sh {params.ref_path} {input} {params.outdir} {params.bbmap}
    """
    
    
  ###############################################################################
 ########################## META TABLE PRODUCTION ##############################
###############################################################################

if config['seurat_obj']!='NA':

  rule filtering_cleaned:
    input:
      ancient(expand(os.path.join(tempo_out, "cleaned/{plate}_{cell_id}.cleaned.fq.gz"), plate=plates, cell_id=cells))

    conda:
      "seurat2"
      
    params:
      patient=config["patient"],
      path=os.path.join(tempo_out, "cleaned"),
      seurat=config["seurat_obj"],
      lab=config["cluster_labs"],
      thr=config["Thr"],
      out=os.path.join(config["working_dir"],config["patient"]),
      label=config["labeling_type"],
      working_dir=config["working_dir"],
      
    output:
      f1= os.path.join(config["working_dir"], config["patient"], "cleaned/samples.tsv"),
      f2= os.path.join(config["working_dir"], config["patient"], "cleaned/units.tsv"),
      f3= os.path.join(config["working_dir"], config["patient"], "cleaned/full_metadata.tsv"),
      
    shell:
      """
      cat {params.path}/*_fracs.txt > {params.path}/fractions.txt
      Rscript scripts/prepare_tables.R {params.path} {params.seurat} '{params.lab}' {params.thr} {params.out} {params.label} {params.working_dir}
      sleep 5
      if [ -s {params.path}/samples.tsv ] && [ -s {params.path}/units.tsv ] && [ -s {params.path}/full_metadata.tsv ]; then
        rsync -av --partial --progress {params.path} {params.out}
      else
        echo output not present
        exit -1
      fi
      """

else:
  rule filtering_cleaned:
    input:
      ancient(expand(os.path.join(config["tmp"], config["patient"], "cleaned/{plate}_{cell_id}.cleaned.fq.gz"), plate=plates, cell_id=cells))
    
    params:
      patient=config["patient"],
      out_dir=os.path.join(config["working_dir"], config["patient"], "cleaned"),
      
    shell:
      '''
      cp {input} -t {params.out_dir}
      '''
   