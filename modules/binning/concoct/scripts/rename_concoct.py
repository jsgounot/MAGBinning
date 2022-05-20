# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-20 13:48:18
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-10-20 13:55:20

import glob, os

fnames = snakemake.params[0]
fnames = glob.glob(fnames)
sample = snakemake.wildcards.sample

for fname in fnames:
	outfile = os.path.join(os.path.dirname(fname), sample + '.bin.' + os.path.basename(fname))
	print (fname, outfile)
	exit()