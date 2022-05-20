# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 16:51:41
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-20 18:14:58

import os
import pandas as pd

rmext = lambda fname : os.path.splitext(fname)[0]
dname = os.path.dirname
bname = os.path.basename

checkm   = snakemake.input.checkm
barrnap  = snakemake.input.barrnap
trnascan = snakemake.input.trnascan
gunc     = snakemake.input.gunc

def isfileempty(fname):
	return os.stat(fname).st_size == 0

def load_tsv(fname, ** kwargs) :
	df = pd.read_csv(fname, sep="\t", ** kwargs)
	df['sample'] = bname(dname(dname(dname(dname(fname)))))
	return df

checkm   = pd.concat((load_tsv(fname, index_col=0) for fname in checkm if not isfileempty(fname)))
barrnap  = pd.concat((load_tsv(fname, index_col=0) for fname in barrnap if not isfileempty(fname)))
trnascan = pd.concat((load_tsv(fname, index_col=0) for fname in trnascan if not isfileempty(fname)))
gunc     = pd.concat((load_tsv(fname, index_col=0) for fname in gunc if not isfileempty(fname)))

trnascan = trnascan[["ID", "tRNAKind", "sample", "tRNAUniqueStrict"]]

# We select the tRNA scan results which display the highest number of tRNA
# this is not the best but the alternative would be to identify
# whether each bin is archea or bacteria, which would need phylogenetic placement
trnascan = trnascan[trnascan.groupby(["sample", "ID"])["tRNAUniqueStrict"].transform(max) == trnascan["tRNAUniqueStrict"]]

# In case of duplicates (# archea == # bacterial), we select the bacterial one
# Change nothing, just make more sense to assign bacteria
trnascan = trnascan.sort_values(["sample", "ID", "tRNAKind"])
trnascan = trnascan.drop_duplicates(subset=["sample", "ID"], keep="last")

df = checkm.merge(barrnap, on=["sample", "ID"], how="left")
df = df.merge(trnascan, on=["sample", "ID"], how="left")

gunc = gunc[['sample', 'genome', 'pass.GUNC']]
gunc.columns = ['sample', 'ID', 'pass.GUNC']
df = df.merge(gunc, on=["sample", "ID"], how="left" )

df["RNAPASS"] = (df["5S_rRNA"] > 0) & (df["16S_rRNA"] > 0) & (df["23S_rRNA"] > 0) & (df["tRNAUniqueStrict"] >= 18)

def get_mimag(row) :
	cstat = row["CheckMStatus"]
	tstat = row["RNAPASS"]

	if cstat == "HIGH" :
		if not tstat : cstat = "NEARCOMPLETE"

	return cstat

df["MIMAG"] = df.apply(get_mimag, axis=1)

outfile = snakemake.output.bins
df.to_csv(outfile, sep="\t", index=False)

# MAGs softlink (MEDIUM AND HIGH)
outdir = snakemake.params.magdir
binner = snakemake.wildcards.binner
os.makedirs(outdir, exist_ok=True)
cwd = os.getcwd()

df = df[df["MIMAG"].isin(("MEDIUM", "HIGH"))]

outfile = snakemake.output.mags
df.to_csv(outfile, sep="\t", index=False)

binlist = snakemake.input.binlist
binpath = {}
for fname in binlist:
	with open(fname) as f:
		for line in f:
			line = line.strip()
			binID = os.path.basename(line)
			binpath[binID] = line

for binID in df['ID']:
	fname = binpath[binID]
	outfname = os.path.join(outdir, binID)
	if not os.path.isfile(outfname) : os.symlink(fname, outfname)

