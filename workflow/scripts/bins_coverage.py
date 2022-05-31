# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-05-27 18:59:54
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-31 14:58:12

import os
import pandas as pd

bincontigs  = snakemake.input.bincontigs
bincoverage = snakemake.input.bincoverage
outfile = snakemake.output[0]

# Read contig files
results = []
with open(bincontigs) as f:
	for line in f :
		binid, contigid = line.strip().split(':')
		binid = os.path.basename(binid)
		contigid = contigid[1:]
		results.append({
			'ID': binid, 'contigid': contigid.split()[0]
			})

df = pd.DataFrame(results)

# Merge coverage info
cov = pd.read_csv(bincoverage, sep=' ', names=['contigid', 'cov', 'npos', 'meancov'])
df = df.merge(cov, on='contigid', how='left')

if not df[df['cov'].isnull()].empty:
	raise Exception('A contig has no coverage ?!', df[df['cov'].isnull()])

df = df.groupby('ID')[['cov', 'npos']].sum().reset_index()
df['meancov'] = df['cov'] / df['npos']

df.to_csv(outfile, sep='\t', index=False)