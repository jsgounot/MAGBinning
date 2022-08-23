# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-08-22 18:05:24
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-08-22 20:56:34

def get_aln_type(config, wc):
	sample = wc.sample
	if 'r1' in config[sample]:
		return 'short'
	elif 'nanopore' in config[sample]:
		return 'long'
	else:
		raise Exception(f'Unable to find either r1/r2 or nanopore for sample {sample}')

def get_aln(config, wc, suffix='.bam'):
	aln = get_aln_type(config, wc)
	return f'data/mapping/{aln}/{wc.sample}/{wc.sample}{suffix}'

def get_aln_vamb(config, wc):
	aln = get_aln_type(config, wc)
	return f'binning/vamb_process/raw/preprocessing/mapping/{wc.sample}.{aln}.bam'

def get_jgi_identity(config, wc):
	aln = get_aln_type(config, wc)
	default = {'short': 97, 'long': 85}
	default = default[get_aln_type(config, wc)]
	ide = config[wc.sample].get('jgi_identity', default)
	return ide