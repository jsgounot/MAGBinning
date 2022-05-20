# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-28 17:18:21
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-10-28 17:20:30

import os
from dereplicate import process

checkmfile = snakemake.input[0]
outdir = os.path.dirname(snakemake.output[0])
partial = bool(snakemake.params['partial'])
makesym = bool(snakemake.params['makesym'])

process(checkmfile, outdir, partial, makesym)