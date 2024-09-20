#!/bin/bash
module load BBMap

hs_38="/g/korbel/zafferani/genome_indexed/human_gn/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
mm_38="/g/korbel/zafferani/genome_indexed/mouse_gn/Mus_musculus.GRCm38.dna.primary_assembly.fa"

##### build index ####
bbsplit.sh build=1 ref_x=$hs_38 ref_y=$mm_38 -Xmx32g
