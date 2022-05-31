# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 16:51:41
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-30 17:53:27

import os
import pandas as pd

rmext = lambda fname : os.path.splitext(fname)[0]
dname = os.path.dirname
bname = os.path.basename

bincov   = snakemake.input.bincov
checkm   = snakemake.input.checkm
barrnap  = snakemake.input.barrnap
trnascan = snakemake.input.trnascan
gunc     = snakemake.input.gunc

def isfileempty(fname):
	return os.stat(fname).st_size == 0

def dfname(fname, fold):
	while fold > 0:
		fname = dname(fname)
		fold -= 1
	return fname

def load_tsv(fname, dfold=4, ** kwargs) :
	df = pd.read_csv(fname, sep="\t", ** kwargs)
	df['sample'] = bname(dfname(fname, dfold))
	return df

bincov   = load_tsv(bincov, dfold=2) 
checkm   = load_tsv(checkm, index_col=0)
barrnap  = load_tsv(barrnap, index_col=0)
trnascan = load_tsv(trnascan, index_col=0)
gunc     = load_tsv(gunc, index_col=0)


bincov = bincov[["sample", "ID", "meancov"]]

trnascan = trnascan[["ID", "tRNAKind", "sample", "tRNAUniqueStrict"]]

# We select the tRNA scan results which display the highest number of tRNA
# this is not the best but the alternative would be to identify
# whether each bin is archea or bacteria, which would need phylogenetic placement
trnascan = trnascan[trnascan.groupby(["sample", "ID"])["tRNAUniqueStrict"].transform(max) == trnascan["tRNAUniqueStrict"]]

# In case of duplicates (# archea == # bacterial), we select the bacterial one
# Change nothing, just make more sense to assign bacteria
trnascan = trnascan.sort_values(["sample", "ID", "tRNAKind"])
trnascan = trnascan.drop_duplicates(subset=["sample", "ID"], keep="last")

df = bincov.merge(checkm, on=["sample", "ID"], how="left")
df = df.merge(barrnap, on=["sample", "ID"], how="left")
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

outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t', index=False)