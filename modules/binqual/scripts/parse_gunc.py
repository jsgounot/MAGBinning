# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-05-20 18:07:08
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-20 18:10:55

import os
import pandas as pd

bname = lambda line: os.path.basename(line.strip())
rmext = lambda fname: os.path.splitext(fname)[0]

binlist = snakemake.input.binlist
gunc = snakemake.input.gunc

df = pd.read_csv(gunc, sep='\t')

with open(binlist) as f:
	mapper = {rmext(bname(fname)): bname(fname) for fname in f}

df['genome'] = df['genome'].apply(lambda genome: mapper[genome])

outfile = snakemake.output[0]
df.to_csv(outfile, sep='\t')