# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-27 08:24:38
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-10-28 17:13:12

from refiner import process

binlists = snakemake.input
outdir = snakemake.params.outdir
threads = snakemake.threads

process(binlists, outdir, decompose=True, ncore=threads)