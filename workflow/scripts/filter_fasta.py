# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-05-31 11:36:31
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-31 11:59:08

from Bio import SeqIO

fname = snakemake.input[0]
minsize = snakemake.params.minsize
records = (record for record in SeqIO.parse(fname, 'fasta') if len(record.seq) >= minsize)
outfile = snakemake.output[0]
SeqIO.write(records, outfile, 'fasta')