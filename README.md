## Another Snakemake binning module

Personal [snakemake](https://snakemake.readthedocs.io/en/stable/) module for metagenomic binning. This pipeline is mostly a translation and cleaner version of the [metaLAS](https://github.com/SeanChenHCY/metaLAS) pipeline, itself very inspired by the [metaWRAP](https://github.com/bxlab/metaWRAP) pipeline. The current pipeline also includes an additional quality module. 

### Pros

* Fine tuning and versatility of the snakemake pipeline
	* Each pipeline part is a module and can be changed and or removed
	* Allow implementation of other tools more easily
	* Time consuming part were tuned in the pipeline to make the whole process faster
* Rewriting of the [metaWRAP binning](https://i.imgur.com/JL665Qo.png) module
	* Direct implementation of original python scripts for snakemake
	* Installing the whole metaWrap (1Gb) can be avoided
* Add a bin quality module
	* Final bins quality are automatically assessed based on the [MIMAG](https://www.nature.com/articles/nbt.3893) standard
	* Chimeric bins identification with [GUNC](https://grp-bork.embl-community.io/gunc/index.html)

###  Pipeline description

#### Summary

* Inputs are an assembly and reads used for this assembly
* Bins data in different ways:
	* Classic approaches:
		* [concoct](https://github.com/BinPro/CONCOCT)
		* [maxbin2](https://sourceforge.net/projects/maxbin2/)
		* [metabat2](https://bitbucket.org/berkeleylab/metabat/src)
		* [vamb](https://github.com/RasmussenLab/vamb) (see specific note)
		* [metabinner](https://github.com/ziyewang/MetaBinner) (not finished, see note)
		* refined - similar to [metawrap](https://github.com/bxlab/metaWRAP) (see specific note)
* Reassembly of the binning results
	* Currently only [unicycler](https://github.com/rrwick/Unicycler) is used to reassemble bins
* Bins quality
	* [CheckM](https://ecogenomics.github.io/CheckM/)
	* [Barrnap](https://github.com/tseemann/barrnap)
	* [tRNAScan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)
	* [GUNC](https://grp-bork.embl-community.io/gunc/index.html)

#### Output

The pipeline will generate and evaluate your bins. At the end, you will have a directory with all bins, a directory with all mags (both softlinked) and two tables with informations for both MAGs and bins. Note that MAGs are just a subset of bins with medium and above quality.

#### Bins quality

Bins quality are based on these data:
* Completness and contamination provided by checkm
* Number of tRNA provided by tRNAScan-SE
* Number of rRNA provided by barnapp
* Chimeric information provided by GUNC
* Bins coverage based on mapping results against the whole 

| MIMAG group   | Completeness | Contamination | # distinct tRNA | # 5S, 16S and 23S rRNA |
| ------------- | ------------ | ------------- | --------------- | ---------------------- |
| LOW           | < 50         | > 10          | -               | -                      |
| MEDIUM        | > 50         | < 10          | -               | -                      |
| NEAR COMPLETE | > 90         | < 5           | -               | -                      |
| HIGH          | > 90         | < 5           | 18              | At least one each      |

Couple of notes:

* `GUNC` result are not used for MIMAG rank. **I highly recommend anyone to either remove or at least carefully examine chimeric MAGs**
* Some studies do not use rank but make analysis base on a score: `completness - 5 x contamination`. MAGs can be selected using a threshold of score > 50.

### Before runnning the pipeline

Some software cannot be completly setup automaticaly:
* [**MANDATORY**] You need to download [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [conda](https://docs.conda.io/en/latest/). [Mamba](https://github.com/mamba-org/mamba) is nice too.
* [**MANDATORY**] You need to download [gunc database](https://grp-bork.embl-community.io/gunc/installation.html) using gunc and provide the database path into the `binning.snk` file.
* [**OPTIONAL**] If you want to run [opera-ms](https://github.com/CSB5/OPERA-MS) assembly, you will need to install the software yourself [using a conda environment](https://github.com/CSB5/OPERA-MS/issues/48) and provide OPERA-MS executable in the `assembly/Snakefile` file.

### The assembly part

A small pipeline is available with this package to generate assemblies which can be run this way:

```
snakemake -s assemblies.snk --configfile config_assembly.json -d res_assembly --use-conda --cores 16 
```

Where `config_assembly.json` is a json configuration file with reads data for each sample looking like this:
```
{
    "sample": {
        "assembler": "megahit",
        "r1": "path.to.read1.fq.gz",
        "r2": "path.to.read2.fq.gz"
    }
}
```

Available assemblers are: `megahit`, `flye`, `canu`, `operams` (more can be easily added). Reads must be gziped before. Once the process finish, you should have a new configuration file at `res_assembly/config_binning.json` which can be used for the snakemake binner pipeline.

### Binning: The quick way

I highly recommend to have a look at the pipeline configuration process and the resources requirement first. Samples informations (short-reads) and assembly must be provided into a configuration file (named `res_assembly/config_binning.json`).

Once you generated a configuration file (use the template), test the pipeline with a dry run:
`snakemake -s binning.snk --configfile res_assembly/config_binning.json -d res_binning --use-conda --cores 18 -np`

You can then run the pipeline this way:
`snakemake -s binning.snk --configfile res_assembly/config_binning.json -d res_binning --use-conda --cores 18 --resources mem_mb=60000 --rerun-incomplete --keep-going`

This will run the pipeline with 18 cores and 60Gb of memory. Note that you can also run snakemake on a cluster, [see the documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

### Pipeline configuration

The usage of each component is automatically determined based on the required output of the pipeline, defined in the main Snakefile `all` rule. It is encoded this way: `binning/{binner}/quality_data/{postbinning}/tables/bins_info.tsv`, where `binnner` can be  one or several elements of the list `('metabat2', 'maxbin2', 'concoct', 'refined')` and `postbinning` can be from `('raw', 'unicycler')`. For example, if you want to use `metabat2` with `unicyler` reassembly, you should write `binning/metabat2/quality_data/unicycler/tables/bins_info.tsv`. Memory and cores definition are also currently hardcoded based on my current experience to what fits the best to me, but you might want to change that for your need (see next section).

### Specific notes

#### Refined binner

The snakemake implements a `refined` binner which mimic [metawrap](https://github.com/bxlab/metaWRAP) pipeline:

* Generate all potential combinations between metabat2, maxbin2, concoct approaches
* Identified the bins showing the best contamination and completeness values
* Dereplicate from these bins potential duplicated contigs

**Important note**: This approach favors low contamination against high completeness. The first step is to merge the **intersection** of contigs from two bins with similar contigs, and look if the intersection is better than the original two bins. The second step is to keep the bin which shows the better contamination/completeness ratio and remove similar bins (based on 80% similarity). This means **you can have contigs overlap between bins** with this method. **On the resource side**, comparison of all possible combinations is based on checkM, which is slow and require some memory (35Gb), making the whole process much slower. Finally, while this script should make the whole process much faster, some thresholds were removed from the original design and you might have more low quality bins with this pipeline (that you can filter later if you want).

#### VAMB

An integration of [VAMB](https://github.com/RasmussenLab/vamb) is available with this pipeline. Note that VAMB is supposed to be run with GPU but current implementation does not include this aspect, which make the process slower but still doable. For large dataset, you might want to run it with GPU. In this case, I recommend launching first snakemake as usual with option `--until vamb_paste_abundances`, this will lead to the creation of all input files requested by VAMB (which need a mapping for each sample). Once done, modify VAMB rule in `workflow/vamb.snk` by adding the `--cuda` option in the `vamb_params` variable. Relaunch the pipeline with option `--until vamb_vamb` and some GPUs. Once finish, you can run one last time the snakemake pipeline as usual.

VAMB does not offer defined post and pre-processing of data, [which differs from usual binners](https://github.com/RasmussenLab/vamb/blob/master/doc/CAMI2.md). In this pipeline, we use arbitrary filters, encoded in the VAMB options with the `vamb_params`: `-o C -m 1000 --minfasta 500000`; and when we filter the output file with a minimum binsize of 100,000bp with the `minsize` parameter. You might want to modify these.

### Notes on the computer resources

The pipeline requires several tools, each having variable needs. Note that it's nearly impossible to define required computer resources since it highly depends of your input data (biological sample, assembly quality and sequencing depth). I record here tools which might be a bottleneck. Overall I recommand using a system with a minimum of 24CPU and at least ~ 100Gb of RAMs (minimum 40Gb).

#### Overall: Don't trust the progress values and add a bit more

The pipeline uses several times [checkpoint](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) to finely adapts its process to generated bins. Therefore, the number of remaining rules will most likely increase during the pipeline execution. <u>I recommend to add a bit more CPU than what you actually have (5/10%)</u> since some tools steps under-use resources which are given to them.

#### The issue with CheckM

[CheckM lineage workflow](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow) uses pplacer to identify the lineage of the bin into a taxonomic tree. This process is memory intensive and scale linearly. In this pipeline, we assume that you might have a lot of bins, coming at different time. This way, we run Checkm individually for each sample. However, this might be suboptimal for low bins quantity and low memory setup, I hope this will not be a big issue for you.

Additionally, I noticed that Checkm crashes on some cluster configuration.

#### Mapping

All binning tools uses (but not only) coverage information which need assembly mapping. This is done using a classical `BWA` or `minimap2` mapping.

#### Binning refinement

The binning refinment tool use `checkM` completness and contamination values on all binning result **and** binning combinations. This means that for a single sample with a refinement from 3 binners, you have 3 + 4 ``checkM`` analyses to do. `checkM` internally uses `pplacer` which is known to be slow and requires [at least 35Gb of memory](https://github-wiki-see.page/m/Ecogenomics/CheckM/wiki/Installation). This can be improved by using additional `pplacer` threads, but the memory requirement is linear to the number of threads involved, so this is not an ideal solutions. If memory is absolutely not an issue for you, you might want to change `checkM` calls directly in the snakemake pipeline at `workflow/binqual.snk`.

#### Reassembly

This step uses [unicycler](https://github.com/rrwick/Unicycler) for each bin and can be particularly computer intensive, especially if you have high depth. This process involves 32 CPUs and `unicycler` calls several softwares, including [spades](https://www.biostars.org/p/267228/) which is known to take a lot of memory for some well covered samples, [especially when it is too multi-threaded](https://www.biostars.org/p/267228/). You might need to reduce the rule to 16 CPUs. I do not recommand using the reassembly mode unless you have a low number of samples (< 20) or unlimited resources.

## Limiting factors and todo list

* Allow long-read mapping when short-reads are not available
* Unicycler sometimes crashes 
* Add more assemblers
