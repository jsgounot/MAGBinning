# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-27 08:24:38
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-07-15 10:55:50

from refiner import process

binlists = snakemake.input
outdir = snakemake.params.outdir
threads = snakemake.threads

# for debug
# blist = ' '.join(binlists)
# cmdline = f'python ../workflow/scripts/refiner.py -binlists {blist} -outdir {outdir} -decompose -threads {threads}'
# print (cmdline)


process(binlists, outdir, decompose=True, ncore=threads)

