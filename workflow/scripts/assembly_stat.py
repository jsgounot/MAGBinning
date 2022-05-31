# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-05-26 09:07:29
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-26 13:18:34

import os
from collections import Counter
from Bio import SeqIO
import numpy as np
import pandas as pd

def nl50(lengths):
	lengths.sort()
	lengths = lengths[::-1]
	
	sumlength = lengths.sum() / 2
	cumsum = lengths.cumsum()

	l50 = np.argmax(cumsum >= sumlength)
	n50 = lengths[l50]
	l50 += 1 # real index

	return n50, l50

def main(fnames):
	
	rdata = []

	for fname in fnames:

		sample = os.path.splitext(os.path.basename(fname))[0]
		assembler = os.path.basename(os.path.dirname(fname))

		fdata = list(SeqIO.parse(fname, 'fasta'))
		lengths = np.array([len(record.seq) for record in fdata])

		# Done in most analysis, including metaquast and most binner
		# they do not consider any sequence lower than 500bp
		flengths = lengths[lengths >= 500]

		base_n50, base_l50 = nl50(lengths)
		flengths_n50, flengths_l50 = nl50(flengths)

		iterseq = [record.seq.upper() for record in fdata]
		gc = sum(seq.count('G') + seq.count('C') for seq in iterseq)
		n = sum(seq.count('N') for seq in iterseq)

		rdata.append({
			'sample': sample, 'assembler': assembler,
			'nseq': len(lengths), 'size': lengths.sum(), 'n50': base_n50, 'l50': base_l50, '#GC': gc, '#N': n, 
			'min500_nseq': len(flengths), 'min500_size': flengths.sum(), 'min500_n50': flengths_n50, 'min500_l50': flengths_l50
			})

	return pd.DataFrame(rdata)

fnames = snakemake.input
outfile = snakemake.output[0]
df = main(fnames)
df.to_csv(outfile, sep='\t')