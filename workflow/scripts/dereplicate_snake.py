# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-28 17:18:21
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-07-18 16:18:45

import os
from dereplicate import process

checkmfile = snakemake.input[0]
outdir = os.path.dirname(snakemake.output[0])
basename = snakemake.params['basename']
partial = bool(snakemake.params['partial'])
makesym = bool(snakemake.params['makesym'])

process(checkmfile, outdir, basename, partial, makesym)