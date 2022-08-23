## Another Snakemake assembly/binning module

Personal [snakemake](https://snakemake.readthedocs.io/en/stable/) module for metagenomic assembly and binning. This pipeline was at first a translation and cleaner version of the [metaLAS](https://github.com/SeanChenHCY/metaLAS) pipeline, itself very inspired by the [metaWRAP](https://github.com/bxlab/metaWRAP) pipeline. Current pipeline is now much more complete. including an optional assembly step, a complete quality module, a reassembly options, additional binners and quality assessment tools.

### Pros

* Flexibility and reliability
  * Include multiple binners and assembly tools
  * Tuned and tested with multiple datasets
  * Deal with both short and long reads data

* A bin quality module
  * Quality of final bins is automatically assessed based on the [MIMAG](https://www.nature.com/articles/nbt.3893) standard
  * Chimeric bins identification with [GUNC](https://grp-bork.embl-community.io/gunc/index.html)
  * Bins coverage is defined by mapping results
  * Everything is saved into a single summary table
* Fine-tuning and versatility of the snakemake pipeline
  * Each pipeline part is a module and can be changed and or removed
  * Time-consuming parts were tuned in the pipeline to make the whole process faster
  * Allow implementation of new tools easily
* Rewriting of the [metaWRAP binning](https://i.imgur.com/JL665Qo.png) module
  * Direct and faster implementation of original python scripts for snakemake
  * Installing the whole metaWrap (1Gb) can be avoided

###  Pipeline description

#### Summary

* Bins data in different ways
  * [concoct](https://github.com/BinPro/CONCOCT)
  * [maxbin2](https://sourceforge.net/projects/maxbin2/)
  * [metabat2](https://bitbucket.org/berkeleylab/metabat/src)
  * [vamb](https://github.com/RasmussenLab/vamb) ([see specific note](https://github.com/jsgounot/MAGBinning#vamb))
  * [metabinner](https://github.com/ziyewang/MetaBinner) (not finished)
  * refined - similar to [metawrap](https://github.com/bxlab/metaWRAP) ([see specific note](https://github.com/jsgounot/MAGBinning#binning-refinement))
* Reassembly of the binning results
  * Currently, only [unicycler](https://github.com/rrwick/Unicycler) can be used to reassemble bins
* Bins quality
  * [CheckM](https://ecogenomics.github.io/CheckM/)
  * [Checkm2](https://github.com/chklovski/CheckM2)
  * [Barrnap](https://github.com/tseemann/barrnap)
  * [tRNAScan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)
  * [GUNC](https://grp-bork.embl-community.io/gunc/index.html)
* Additional informations
  * MIMAG status 
  * Bins coverage


#### Bins quality

Bins quality is based on these data:
* Completeness and contamination provided by checkm(2)
* Number of tRNA provided by tRNAScan-SE
* Number of rRNA provided by barnapp
* Chimeric information provided by GUNC
* Bins coverage based on mapping results against <u>all assembly contigs</u>

| MIMAG group   | Completeness | Contamination | # distinct tRNA | # 5S, 16S and 23S rRNA |
| ------------- | ------------ | ------------- | --------------- | ---------------------- |
| LOW           | < 50         | > 10          | -               | -                      |
| MEDIUM        | > 50         | < 10          | -               | -                      |
| NEAR COMPLETE | > 90         | < 5           | -               | -                      |
| HIGH          | > 90         | < 5           | 18              | At least one each      |

Couple of notes:

* `GUNC` results are not used for MIMAG rank. **I highly recommend anyone to either remove or at least carefully examine chimeric MAGs**
* Some studies do not use rank but make analysis base on a score: `completness - 5 x contamination`. MAGs can be selected using a threshold of score > 50.

### Setup, conda environments and local variables

While it's in theory possible to run the pipeline without conda, using conda is highly recommended. This pipeline use an extensive amount of software with potential conflict between them. 

* [**MANDATORY**] You need to download [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [conda](https://docs.conda.io/en/latest/). [Mamba](https://github.com/mamba-org/mamba) is nice too.
* [**MANDATORY**] You need to download [gunc database](https://grp-bork.embl-community.io/gunc/installation.html) using gunc and provide the database path into the main snakemake file (see `megahit_metabat2.snk`).
* [**MANDATORY**] If you plan to use `checkm2`, you need to install it yourself and update its path in the main snakemake file (see `megahit_metabat2.snk`).
* [**OPTIONAL**] If you want to run [opera-ms](https://github.com/CSB5/OPERA-MS) assembly, you will need to install the software yourself [using a conda environment](https://github.com/CSB5/OPERA-MS/issues/48) and provide OPERA-MS executable in the `assembly/Snakefile` file.

By default, I set up 7 conda environments (which can be found in `workflow/envs`). The conda environment to use for each software is provided within the file `condaenvs.json` . By default, snakemake will automatically downloads and sets up these environment, but you can change on environment with one of your by editing the json file.

### Binning: The quick way

I highly recommend having a look at the pipeline configuration process and the resources requirement first. Sample information (short-reads) and assembly must be provided in a configuration file (named `config_binning.json`).

Once you generated a configuration file (use the template), test the pipeline with a dry run:
`snakemake -s binning.snk --configfile config_binning.json -d res_binning --use-conda --cores 18 -np`

You can then run the pipeline this way:
`snakemake -s binning.snk --configfile res_assembly/config_binning.json -d res_binning --use-conda --cores 18 --resources mem_mb=60000 --rerun-incomplete --keep-going`

This will run the pipeline with 18 cores and 60Gb of memory.

### With an assembly

Most of the assemblers are not completely set up in the workflow, this part needs to be finished. While the default `binning.snk` workflow expects to find an assembly in the configuration file, you can easily edit this file to include an assembler as provided with `megahit_metabat2.snk` example. Here is the important line in the beginning of the script:

```python
for sample, sdata in config.items():
    sdata['assembly'] = f'assembly/completed/megahit/{sample}.fa'
```

### Pipeline configuration

The usage of each component is automatically determined based on the required output of the pipeline, defined in the main Snakefile `all` rule. It is encoded this way: `binning/{binner}/quality_data/{postbinning}/tables/bins_info.tsv`, where `binnner` can be  one or several elements of the list `('metabat2', 'maxbin2', 'concoct', 'refined')` and `postbinning` can be from `('raw', 'unicycler')`. For example, if you want to use `metabat2` with `unicyler` reassembly, you should write `binning/metabat2/quality_data/unicycler/tables/bins_info.tsv`. Memory and cores amount is also currently hardcoded based on my current experience to what fits the best to me, but you might want to change that for your need (see next section).

### Working with nanopore data

Binning algorithm can deal with both short-reads and nanopore long-reads data, just use `nanopore` as key instead of `r1` and `r2` in your configuration file. You might as well redefined the identity threshold used for the coverage extraction by `jgi_summary` ([see here](https://bitbucket.org/berkeleylab/metabat/issues/142/binning-with-nanopore-data)). This threshold is set to **85** for nanopore data while I keep the script default identity of **97** for short-read data. You can update this value per sample using the key `jgi_identity` in your configuration file as such:

```json
{
	"MYSAMPLE": {
		"assembly": "/path/to/assembly.fa",
		"nanopore": "/path/to/longread.fq.gz",
		"jgi_identity": 90
	}
}
```

If you have both short and long-read, short-reads are prioritized and long-reads are ignored.

### Specific notes

#### Checkm or checkm2

[`Checkm2`](https://github.com/chklovski/CheckM2) is available and can be used flawlessly as a replacement of `checkm`. Note that the two softwares does not produce exactly the same results based on my tests. One major technical advantage of `checkm2` over `checkm` is its independence over `pplacer`, which make the former much faster and memory efficient. This might be especially useful when using the `refined` binner (see below). To use `checkm2` you need to install yourself the conda environment, link the binary in the main snakefile and setup the flag in the internal configuration dictionary (see `megahit_metabat2.snk` for an example).

#### Refined binner

The snakemake implements a `refined` binner that mimic [metawrap](https://github.com/bxlab/metaWRAP) pipeline:

* Generate all potential combinations between metabat2, maxbin2, concoct approaches
* Identified the bins showing the best contamination and completeness values
* Dereplicate from these bins potential duplicated contigs

**Important note**: This approach favors low contamination against high completeness. The first step is to merge the **intersection** of contigs from two bins with similar contigs, and look if the intersection is better than the original two bins. The second step is to keep the bin that shows the better contamination/completeness ratio and remove similar bins (based on 80% similarity). This means **you can have contigs overlap between bins** with this method. **On the resource side**, the comparison of all possible combinations is based on checkM, which is slow and requires some memory (35Gb), making the whole process much slower. Finally, while this script should make the whole process much faster, some thresholds were removed from the original design and you might have more low-quality bins with this pipeline (that you can filter later if you want).

#### VAMB

An integration of [VAMB](https://github.com/RasmussenLab/vamb) is available with this pipeline. Note that VAMB is supposed to be run with GPU but the current implementation does not include this aspect, which makes the process slower but still doable. For a large dataset, you might want to run it with GPU. In this case, I recommend launching first snakemake as usual with option `--until vamb_paste_abundances`, this will lead to the creation of all input files requested by VAMB (which need a mapping for each sample). Once done, modify VAMB rule in `workflow/vamb.snk` by adding the `--cuda` option in the `vamb_params` variable. Relaunch the pipeline with option `--until vamb_vamb` and some GPUs. Once finished, you can run one last time the snakemake pipeline as usual.

VAMB does not offer defined post and pre-processing of data, [which differs from usual binners](https://github.com/RasmussenLab/vamb/blob/master/doc/CAMI2.md). In this pipeline, we use arbitrary filters, encoded in the VAMB options with the `vamb_params`: `-o C -m 1000 --minfasta 500000`; and when we filter the output file with a minimum binsize of 100,000bp with the `minsize` parameter. You might want to modify these.

### Notes on the computer resources

The pipeline requires several tools, each having variable needs. Note that it's nearly impossible to define required computer resources since it highly depends on your input data (biological sample, assembly quality and sequencing depth). I record here tools that might be a bottleneck. Overall I recommend using a system with a minimum of 24CPU and at least ~ 100Gb of RAMs (**minimum 40Gb**). <u>I recommend adding a bit more CPU than what you  have (5/10%)</u> since some tools steps under-use resources which are given to them.

#### Overall: Don't trust the progress values

The pipeline uses several times [checkpoint](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) to finely adapts its process to generated bins. Therefore, the number of remaining rules will most likely increase during the pipeline execution. 

#### The issue with CheckM

[CheckM lineage workflow](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow) uses pplacer to identify the lineage of the bin into a taxonomic tree. This process is memory intensive and scales linearly. In this pipeline, we assume that you might have a lot of bins, coming at different times. This way, we run Checkm individually for each sample. However, this might be suboptimal for low bins quantity and low memory setup, I hope this will not be a big issue for you. I noticed that Checkm crashes on some cluster configurations.

#### Bins refinement - A memory expensive process

The binning refinement tool uses `checkM` completeness and contamination values on all binning result **and** binning combinations. This means that for a single sample with refinement from 3 binners, you have 3 + 4 ``checkM`` analyses to do. `checkM` internally uses `pplacer` which is known to be slow and requires [at least 35Gb of memory](https://github-wiki-see.page/m/Ecogenomics/CheckM/wiki/Installation). This can be improved by using additional `pplacer` threads, but the memory requirement is linear to the number of threads involved, so this is not an ideal solution. If memory is absolutely not an issue for you, you might want to change `checkM` calls directly in the snakemake pipeline at `workflow/binqual.snk`.

#### Reassembly - A CPU expensive process

This step uses [unicycler](https://github.com/rrwick/Unicycler) for each bin and can be particularly intensive, especially if you have high depth. This process involves 32 CPUs and `unicycler` calls several software, including [spades](https://www.biostars.org/p/267228/) which is known to take a lot of memory for some well-covered samples, [especially when it is too multi-threaded](https://www.biostars.org/p/267228/). You might need to reduce the rule to 16 CPUs. I do not recommend using the reassembly mode unless you have a low number of samples (< 20) or unlimited resources.

#### BUG: Out of jobs ready to be started, but not all files built yet

You might hit this message, linked [to this issue](https://github.com/snakemake/snakemake/issues/823), usually almost at the end of the pipeline. This is not a big issue, and you can just run the pipeline again with the same arguments, snakemake will continue where he stopped and finish the last rules.

## Limiting factors and todo list

* Allow long-read mapping when short-reads are not available
* Unicycler sometimes crashes 
* Add more assemblers
