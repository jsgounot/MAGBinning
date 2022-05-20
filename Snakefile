import os, glob

configfile: "config.json"

BINNERS = ('metabat2', 'maxbin2', 'concoct')

rule all:
    input:
        #expand("binning/{binner}/table/bins_info.tsv", binner=BINNERS),
        #expand("binning/{binner}/table/mags_info.tsv", binner=BINNERS)
        #'binning/refined/FMT009_3/dereplicate/checkm/results.tsv'
        #'binning/refined/quality_data/unicycler/tables/mags_info.tsv'
        #directory('binning/refined/FMT005_1/unicycler/binlist.txt'),
        #'binning/maxbin2/quality_data/raw/tables/bins_info.tsv',
        'binning/maxbin2/quality_data/unicycler/tables/bins_info.tsv',
        #'binning/refined/quality_data/raw/tables/bins_info.tsv',
        #'binning/refined/quality_data/unicycler/tables/bins_info.tsv'

# Import modules

module mapping:
    snakefile: "modules/mapping/Snakefile"
    config: config

use rule * from mapping as mapping_*

module concoct:
    snakefile: "modules/binning/concoct/Snakefile"

use rule * from concoct as concoct_*

module maxbin2:
    snakefile: "modules/binning/maxbin2/Snakefile"

use rule * from maxbin2 as maxbin2_*

module metabat2:
    snakefile: "modules/binning/metabat2/Snakefile"

use rule * from metabat2 as metabat2_*

module binrefiner:
    snakefile: "modules/binrefiner/Snakefile"
    config: config

use rule * from binrefiner as binrefiner_*

module binqual:
    snakefile: "modules/binqual/Snakefile"
    config: config

use rule * from binqual as binqual_*

module reads_retrieval:
    snakefile: "modules/postbinning/reads_retrieval/Snakefile"
    config: config

use rule * from reads_retrieval as reads_retrieval_*

module unicycler:
    snakefile: "modules/postbinning/unicycler/Snakefile"
    config: config

use rule * from unicycler as unicycler_*