# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-08-22 18:05:24
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-10-20 17:09:26

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

def get_checkm_dir(config):
    checkm_version = str(config['softparams']['softdb']['checkm_version']).strip()
    checkm_dir = {'1': 'checkm', '2': 'checkm2'}.get(checkm_version, None)
    if checkm_dir is None: raise Exception('Checkm version must be 1 or 2')
    return checkm_dir