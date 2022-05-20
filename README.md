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
		* [CheckM](https://ecogenomics.github.io/CheckM/)
		* [Barrnap](https://github.com/tseemann/barrnap)
		* [tRNAScan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)

###  Summary of the pipeline

* Inputs are an assembly and reads used for this assembly

* Bins data in different ways:
	* Classic approaches:
		* [concoct](https://github.com/BinPro/CONCOCT)
		* [maxbin2](https://sourceforge.net/projects/maxbin2/)
		* [metabat2](https://bitbucket.org/berkeleylab/metabat/src)
	* Refined module:
		* Generate all potential combinations between classic approaches
		* Identified the bins showing the best contamination and completness values
		* Dereplicate from these bins potential duplicated contigs
* Reassembly of the binning results
	* Currently only [unicycler](https://github.com/rrwick/Unicycler) is used to reassemble bins
* Bins quality
		* [CheckM](https://ecogenomics.github.io/CheckM/)
		* [Barrnap](https://github.com/tseemann/barrnap)
		* [tRNAScan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)

### The quick way

I highly recommend to have a look at the pipeline configuration process and the resources requirement first.

Once you generated a configuration file (use the template), test the pipeline with a dry run:
`snakemake -np`

If no errors are shown, run the pipeline like this:
`snakemake --cores 24 --use-conda --resources mem_mb=60000`

This will run the pipeline with 24 cores and 60Gb of memory.

You might want to add these options:
`snakemake --cores 24 --use-conda --resources mem_mb=60000 --rerun-incomplete --keep-going`

Note that you can also run snakemake on a cluster, [see the documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

### Pipeline configuration

The usage of each component is automatically determined based on the required output of the pipeline, defined in the main Snakefile `all` rule. It is encoded this way: `binning/{binner}/quality_data/{postbinning}/tables/bins_info.tsv`, where `binnner` can be  one or several elements of the list `('metabat2', 'maxbin2', 'concoct', 'refined')` and `postbinning` can be from `('raw', 'unicycler')`. For example, if you want to use `metabat2` with `unicyler` reassembly, you should write `binning/metabat2/quality_data/unicycler/tables/bins_info.tsv`. Memory and cores definition are also currently hardcoded based on my current experience to what fits the best to me, but you might want to change that for your need (see next section).

### Notes on the computer resources

The pipeline requires several tools, each having variable needs. Note that it's nearly impossible to define required computer resources since it highly depends of your input data (biological sample, assembly quality and sequencing depth). I record here tools which might be a buttleneck. Overall I recommand using a system with a minimum of 16CPU and at least ~ 100Gb of RAMs.

#### Overall: Don't trust the progress values and add a bit more

The pipeline uses several times [checkpoint](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution)  to finely adapted its process to formed bins. Therefore, the number of remaining rules will most likely increase during the pipeline execution. 

I also recommend to add a bit more CPU than what you actually have (5/10%) since some tools steps under-use resources which are given to them.

#### Mapping

All binning tools uses (but not only) coverage information which need assembly mapping. This is done using a classical `BWA` mapping which can take some time and CPU power (default 8 CPU).

#### Binning refinement

The binning refinment tool use `checkM` completness and contamination values on all binning result **and** binning combinations. This means that for a single sample with a refinement from 3 binners, you have 3 + 4 ``checkM`` analyses to do. `checkM` is known to be slow and requires [at least 35Gb of memory](https://github-wiki-see.page/m/Ecogenomics/CheckM/wiki/Installation). This can be improved by using additional `pplacer` threads, but the memory requirement is linear to the number of threads involved, so this is not an ideal solutions. If memory is absolutely not an issue for you, you might want to change `checkM` calls directly in the scripts.

#### Reassembly

This step is using [unicycler](https://github.com/rrwick/Unicycler)  for each bin individually and can be particularly computer intensive, especially if you have high depth. This step involved 32 CPUs and `unicycler` calls several softwares, including [spades](https://www.biostars.org/p/267228/) which is known to take a lot of memory for some well covered samples, [especially when it is too multi-threaded](https://www.biostars.org/p/267228/). You might need to reduce the rule to 16 CPUs. 

#### Quality

Involved `checkM` (see bin refinement), but just one call per sample. While `barrnap` is called once per sample, `tRNA-Scan` is launched individually for each bin.


