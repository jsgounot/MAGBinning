# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-28 17:13:04
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-10-28 17:30:51

import os
from consolidate import process

inputs = list(snakemake.input)
outdir = os.path.dirname(snakemake.output[0])
overlap = snakemake.params['overlap']

process(inputs, outdir, overlap)