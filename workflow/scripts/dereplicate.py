# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-28 15:24:07
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-31 13:33:08

"""
Adapted version of the metaWrap script: https://github.com/bxlab/metaWRAP/blob/master/bin/metawrap-scripts/dereplicate_contigs_in_bins.py
"""

import argparse, os, pathlib, shutil
import pandas as pd
from Bio import SeqIO

def load_contigs(df):
    for binid, fname in zip(df['NewID'], df['binfname']):
        for record in SeqIO.parse(fname, 'fasta'):
            yield binid, record.id

def make_abspath_check(fname):
    fname = pathlib.Path(fname)
    fname = str(fname.absolute())

    if not os.path.isfile(fname):
        raise OSError('Unable tor retrieve fname: ' + fname)

    return fname

def touch(fname):
    with open(fname, 'w') as f:
        pass

def process(checkmfile, outdir, partial=True, makesym=True):
    df = pd.read_csv(checkmfile, sep='\t')

    if df.empty:
        outfile = os.path.join(outdir, 'binlist.txt')
        touch(outfile)
        exit(0)

    df['score'] = df['Completness'] - 5 * df['Contamination'] + df['Contig_N50'] * 1e-10
    
    directory = os.path.dirname(checkmfile)
    df['binfname'] = df['NewID'].apply(lambda binid: os.path.join(directory, 'bin.' + str(binid) + '.fa'))

    sdf = pd.DataFrame([{'NewID': binid, 'contig': contig} 
        for binid, contig in load_contigs(df)])

    count = sdf.groupby('contig').size().to_dict()
    sdf['count'] = sdf['contig'].map(count)

    sdf = sdf[sdf['count'] > 1]
    sdf = sdf.merge(df[['NewID', 'score']], on='NewID', how='left')

    if partial:
        # We remove duplicated contigs from bins which display
        # the lowest score. To do this we sort bins by score (from the highest to the lowest score),
        # extract the best bin for each contig by dropping duplicates and keeping last (best score)
        # and use this info to retrieve all other rows (bad scores)
        sdf = sdf.sort_values(['contig', 'score'])
        sdf = sdf[~ sdf.index.isin(sdf.drop_duplicates('contig', keep='last').index)]

    # if not partial, we keep the sdf which contains all the duplicates (count > 1)
    # we can now remove data
    torm = sdf.groupby('NewID')['contig'].apply(set)
    makefname = lambda binid: os.path.join(outdir, 'bin.' + str(binid) + '.fa')
    outfiles = []

    for binid, fname in zip(df['NewID'], df['binfname']):
        outfname = makefname(binid)

        if binid in torm:
            fdata = SeqIO.parse(fname, 'fasta')
            fdata = list(record for record in fdata if record.id not in torm[binid])
            if not fdata: continue
            SeqIO.write(fdata, outfname, 'fasta')

        elif makesym:
            fname = make_abspath_check(fname)
            os.symlink(fname, outfname)

        else:
            shutil.copyfile(fname, outfname)

        outfiles.append(make_abspath_check(outfname))

    outfile = os.path.join(outdir, 'removed.contigs.tsv')
    sdf.to_csv(outfile, sep='\t', index=False)

    outfile = os.path.join(outdir, 'binlist.txt')
    with open(outfile, 'w') as f:
        f.write('\n'.join(outfiles) + '\n')

# -------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-checkmfile',
                        required=True,
                        help='Path to checkmfile result')

    parser.add_argument('-outdir',
                        required=True,
                        help='Path to outdir')

    parser.add_argument('-partial',
                        action='store_true',
                        help='Complete dereplication, will remove duplicated contigs in all bins. Default partial: Keep contig in the best bin.')

    parser.add_argument('-makesym',
                        action='store_true',
                        help='Will generate symlink instead of copying file for untouched bins')

    args = parser.parse_args()

    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    process(args.checkmfile, outdir, args.partial, args.makesym)