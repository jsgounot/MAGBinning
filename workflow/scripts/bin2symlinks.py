# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-07-27 09:44:54
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-07-27 17:12:57

# This one is just a safe replacement of a bash loop
# I've encounter a couple of issues with bash so I make it easier here

import os

binlist = snakemake.input[0]
outdir  = snakemake.output[0]
ext = snakemake.params.get('ext', None)

os.makedirs(outdir, exist_ok=True)

with open(binlist) as f:
	for fname in f:
		fname = fname.strip()
		if not fname: continue
		if not os.path.isfile(fname):
			raise IOError(f'File not found {fname}')

		outfile = os.path.join(outdir, os.path.basename(fname))

		if ext is not None:
			outfile = os.path.splitext(outfile)[0] + '.' + ext

		if os.path.isfile(outfile): 
			os.remove(outfile)
		
		os.symlink(fname, outfile)