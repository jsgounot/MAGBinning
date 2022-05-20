# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-11-01 14:43:02
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-11-01 21:17:46

import os

binlist = snakemake.input[0]
bins = list(snakemake.input[1:])

output = snakemake.output[0]
outdir = snakemake.output[1]

with open(binlist) as f:
	binlist = list(f)

print (bins)
exit()

with open(output, 'w') as f:
	f.write('\n'.join(inputs))