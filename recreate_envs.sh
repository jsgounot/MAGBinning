# @Author: jsgounot
# @Date:   2022-05-20 15:39:23
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-20 17:38:55

mamba create -n concoct -c bioconda concoct
conda activate concoct
conda env export > modules/binning/concoct/envs/concoct.yaml
conda activate snakemake

mamba create -n maxbin2 -c bioconda maxbin2
conda activate maxbin2
conda env export > modules/binning/maxbin2/envs/maxbin2.yaml
conda activate snakemake

mamba create -n metabat -c bioconda metabat2
conda activate metabat
conda env export > modules/binning/metabat2/env/metabat2.yaml
conda activate snakemake

mamba create -n gunc -c bioconda gunc
conda activate gunc
conda env export > modules/binqual/env/gunc.yaml
conda activate snakemake

mamba create -n checkm -c bioconda checkm-genome
conda activate checkm
conda env export > modules/binqual/env/checkm.yaml
conda env export > modules/binrefiner/env/checkm.yaml
conda activate snakemake

mamba create -n biopy -c conda-forge biopython
conda activate biopy
conda env export > modules/binrefiner/env/biopy.yaml
conda activate snakemake

mamba create -n mapping -c bioconda bwa minimap2 samtools
conda activate mapping
conda env export > modules/mapping/env/mapping.yaml
conda activate snakemake

mamba create -n bamtools -c bioconda picard samtools
conda activate bamtools
conda env export > modules/postbinning/reads_retrieval/env/bamtools.yaml
conda activate snakemake

mamba create -n unicycler -c bioconda unicycler
conda activate unicycler
conda env export > modules/postbinning/unicycler/env/unicycler.yaml
conda activate snakemake