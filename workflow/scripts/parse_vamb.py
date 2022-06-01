# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-05-31 17:58:30
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-06-01 10:36:14

'''
In this script:
1. We check that we have the correct bin association
2. Provide a stat file with bin info
3. Split bins for each strain
'''

import gzip, os, shutil
import pandas as pd
from Bio import SeqIO

cluster = snakemake.input.cluster
contigs = snakemake.input.contigs[0]
indcons = snakemake.input.indcons 
minsize = snakemake.params.minsize
samples = snakemake.params.samples
binsdir = snakemake.params.binsdir
outfile = snakemake.output[0]

extract_original_contig_name = eocn = lambda contig: contig[contig.index('C')+1:]

with gzip.open(contigs, 'rt') as f:
	fdata = list(SeqIO.parse(f, 'fasta'))
	fleng = {record.id: len(record.seq) for record in fdata}

fdataind = {str(fidx + 1): {record.id: len(record.seq)
	for record in SeqIO.parse(fasta, 'fasta')}
	for fidx, fasta in enumerate(indcons)}

df = pd.read_csv(cluster, sep='\t', names=['cluster', 'contig'])
df['main_length'] = df['contig'].map(fleng)

df['sample_index'] = df['contig'].apply(lambda contig: contig[1:contig.index('C')])
df['contig_id'] = df['contig'].apply(eocn)
df['length_id'] = df.apply(lambda row: fdataind[row['sample_index']][row['contig_id']], axis=1)

# very important to check this, otherwise we might associated the wrong
# fasta to one sample. It's not completly bullet proof though
assert df[df['length_id'] != df['main_length']].empty

sdf = df.groupby(['cluster', 'sample_index'])['length_id'].agg(['size', 'sum']).reset_index()
sdf['sample_name'] = sdf['sample_index'].apply(lambda idx: samples[int(idx) - 1])
sdf = sdf[['cluster', 'sample_index', 'sample_name', 'size', 'sum']]
sdf.columns = ['vamb_cluster', 'sample_index', 'sample_name', '# contigs', 'total length']

sdf = sdf[sdf['total length'] >= minsize]
sdf.to_csv(outfile, sep='\t', index=False)

for record in fdata:
	record.vamb_id = record.id
	record.id = eocn(record.id)

for sinfo, ssdf in sdf.groupby(['sample_name', 'sample_index']):
	sample, sindex = sinfo
	sdata = fdataind[sindex]
	
	outdir = os.path.join(binsdir, sample, 'raw/bins')
	if os.path.isdir(outdir): shutil.rmtree(outdir)
	os.makedirs(outdir)

	for cluster in ssdf['vamb_cluster']:
		outfile = os.path.join(outdir, sample + '.bin.' + cluster + '.fa')
		# print (f'Write bin in {outfile}')

		contigs = set(df[(df['cluster'] == cluster) & (df['sample_index'] == sindex)]['contig'])
		sdata = [record for record in fdata if record.vamb_id in contigs]

		assert len(sdata) == len(contigs)
		SeqIO.write(sdata, outfile, 'fasta')
