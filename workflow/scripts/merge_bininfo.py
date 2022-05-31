# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-05-30 17:50:06
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-30 19:13:20

import pandas as pd
import os

def touch(fname):
	with open(fname, 'w') as f:
		f.write('')

bininfos = snakemake.input.bininfos
notempty = lambda fname: os.stat(fname).st_size != 0
df = [pd.read_csv(fname, sep='\t') for fname in bininfos if notempty(fname)]

binstat = snakemake.output.bins
magstat = snakemake.output.mags

if not df:
	# we only have empty files <- No bins at all
	touch(binstat)
	touch(magstat)
	exit(0)

df = pd.concat(df)

# Load bin paths
binlists = snakemake.input.binlists
binpaths = {}
for fname in binlists:
	with open(fname) as f:
		for line in f:
			line = line.strip()
			binID = os.path.basename(line)
			binpaths[binID] = line

def link_bins(binIDs, binpaths, outdir):
	os.makedirs(outdir, exist_ok=True)
	for binID in binIDs:
		fname = binpaths[binID]
		outfname = os.path.join(outdir, binID)
		if not os.path.isfile(outfname) : os.symlink(fname, outfname)

# Bins tab and links
df.to_csv(binstat, sep="\t", index=False)
outdir = snakemake.params.bindir
link_bins(df['ID'], binpaths, outdir)

# MAGs tab and links
df = df[df["MIMAG"].isin(("MEDIUM", "HIGH"))]
df.to_csv(magstat, sep="\t", index=False)
outdir = snakemake.params.magdir
link_bins(df['ID'], binpaths, outdir)