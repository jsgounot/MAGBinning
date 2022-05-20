# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-26 22:05:51
# @Last Modified by:   jsgounot
# @Last Modified time: 2021-10-28 17:39:32

"""
Part of the MetaWRAP binrefiner module: https://github.com/bxlab/metaWRAP
Adapted from this script: https://github.com/bxlab/metaWRAP/blob/master/bin/metawrap-scripts/consolidate_two_sets_of_bins.py

Summary:
This script look for the best bins among a group of bins with an iterative process
1. Will pairwise compare bins and define two bins are similar if they share x (default 80%) similarity (see note A)
2. Will select the bin with the best completness and contamination based on the formula: completness - contamination * 5

Note:
A. Similar to the refiner script, this process assumes that bins are coming from the same assembly and that contig sharing the same
   name are similar.
B. This script take as input a processed checkm file. See parse_checkm.py script for this.

Improvement:
- One call for all combinations of the iteration process. No overload of data
- Panda dataframe comparison instead of nested dictionaries
- No multiple copies and renaming of the data
- Extract and yield a reconstructed checkM datafile for the resulting bins
"""

import os, pathlib
import pandas as pd
from Bio import SeqIO

dname = os.path.dirname

def extract_contigs(checkm_df):
    binsinfo = []
    make_bin_path = lambda binid, checkrpath: os.path.join(dname(checkrpath), 'inputs', binid + '.checkminput')

    for binid, sfname in zip(checkm_df['ID'], checkm_df['Source']):
        binfile = make_bin_path(binid, sfname)
        
        if not os.path.isfile(binfile):
            raise OSError('CheckM input file not found: ' + binfile)

        for record in SeqIO.parse(binfile, 'fasta') :
            binsinfo.append({'binID': binid, 'contig': record.id, 'length': len(record.seq)})

    return pd.DataFrame(binsinfo)

def load_initial_data(checkm_result, contig_data, checkm_data):
    print ('Reading checkM file: ' + checkm_result)

    df = pd.read_csv(checkm_result, sep='\t', index_col=0)
    sdf = df[(df['Contamination'] < 10) & (df['Completness'] > 70)]
    sdf['Source'] = checkm_result

    print ('%i bins total and %i passing the filters' %(df['ID'].nunique(), sdf['ID'].nunique()))
    checkm_data[checkm_result] = sdf
    contig_data[checkm_result] = extract_contigs(sdf)

def extract_best(row):
    names = {'ID': 'binID', 'Completness': 'completness', 'Contamination': 'contamination', 'Source': 'source'}
    suffix = '_B' if row['better_B'] else '_A'
    return pd.Series({name: row[cname + suffix] for name, cname in names.items()})

def pairwise_process(checkm_result_A, checkm_result_B, contig_data, checkm_data, overlap, rnum):
    cMdfA, cIdfA = checkm_data[checkm_result_A], contig_data[checkm_result_A]
    cMdfB, cIdfB = checkm_data[checkm_result_B], contig_data[checkm_result_B]

    mdf = cIdfA.merge(cIdfB, on='contig', how='inner', suffixes=('_A', '_B'))
    mdf = mdf.groupby(['binID_A', 'binID_B'])[['length_A', 'length_B']].sum().reset_index()
        
    # I guess lengthA == lengthB but I keep what have been done previously
    # This would not be the case in we have contig with the same name but different size
    # Which would mean coming from different assembly and / or have been through post-processing

    mdf['binSize_A'] = mdf['binID_A'].map(cIdfA.groupby('binID')['length'].sum())
    mdf['binSize_B'] = mdf['binID_B'].map(cIdfB.groupby('binID')['length'].sum())

    mdf['ratio_A'] = 100 * mdf['length_A'] / mdf['binSize_A']
    mdf['ratio_B'] = 100 * mdf['length_B'] / mdf['binSize_B']
    mdf['bestRatio'] = mdf[['ratio_A', 'ratio_B']].max(axis=1)

    # We prune based on overlap values
    mdf = mdf[mdf['bestRatio'] >= 80]   

    # We keep the bins which are not overlapping to add them later
    cMdfA_notfound = cMdfA[~ cMdfA['ID'].isin(mdf['binID_A'])][['ID', 'Completness', 'Contamination', 'Source']]
    cMdfB_notfound = cMdfB[~ cMdfB['ID'].isin(mdf['binID_B'])][['ID', 'Completness', 'Contamination', 'Source']]

    # Adding checkM values in the table for bins set A
    sdf = cMdfA[['ID', 'Completness', 'Contamination', 'Source']]
    sdf.columns = ['binID_A', 'completness_A', 'contamination_A', 'source_A']
    mdf = mdf.merge(sdf, on='binID_A')
    mdf['score_A'] = mdf['completness_A'] - mdf['contamination_A'] * 5

    # Adding checkM values in the table for bins set B
    sdf = cMdfB[['ID', 'Completness', 'Contamination', 'Source']]
    sdf.columns = ['binID_B', 'completness_B', 'contamination_B', 'source_B']
    mdf = mdf.merge(sdf, on='binID_B')
    mdf['score_B'] = mdf['completness_B'] - mdf['contamination_B'] * 5

    # We select the overlapping bin for each binID_A which have the best score
    # and compare this value to the binA score to select which bin we want to conserve
    mdf = mdf[mdf['score_B'] == mdf.groupby('binID_A')['score_B'].transform('max')]
    mdf['better_B'] = mdf['score_B'] > mdf['score_A']

    # We then produce a new dataframe based on result founds before
    # We rename all best results as bin_A to be used next in the iterative process
    mdf = mdf.apply(extract_best, axis=1)

    # We add bins which pass the checkm threshold but were not include because they
    # were not found in the overlapping process (both sides)
    mdf = pd.concat([mdf, cMdfA_notfound, cMdfB_notfound])
    print ('Round %i: %i bins cherry-picked' %(rnum, len(mdf)))

    # We create a new dataset with the new results
    # This could be faster by just retrieving previous informations (TODO)
    checkm_data['current'] = mdf
    contig_data['current'] = extract_contigs(mdf)

def reconstruct_checkm_results(finaldf, checkm_data):
    found = finaldf.groupby('Source')['ID'].apply(set).to_dict()

    def extract_sub(source, binids, checkm_data):
        sdf = checkm_data[source]
        return sdf[sdf['ID'].isin(binids)]

    sdf = pd.concat((extract_sub(source, binids, checkm_data)
        for source, binids in found.items()))

    fdf = finaldf[['Source', 'ID', 'NewID']]
    sdf = sdf.merge(fdf, on=['Source', 'ID'], how='outer')
    return sdf.sort_values('NewID')

def save_contigs(finaldf, outdir):
    
    def make_bin_path(binid, checkrpath):
        fname = os.path.join(dname(checkrpath), 'inputs', binid + '.checkminput')
        fname = pathlib.Path(fname)
        fname = str(fname.absolute())

        if not os.path.isfile(fname):
            raise OSError('Unable tor retrieve fname: ' + fname)

        return fname

    idx = 1

    for binid, sfname, newid in zip(finaldf['ID'], finaldf['Source'], finaldf['NewID']):
        binfile = make_bin_path(binid, sfname)
        outfile = os.path.join(outdir, 'bin.' + str(newid) + '.fa')
        os.symlink(binfile, outfile)
        idx += 1

def process(fnames, outdir, overlap=80):
    contig_data = {}
    checkm_data = {}

    if not len(fnames) > 1:
        raise ValueError('fnames must be a list with at least two checkM paths')

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    for fname in fnames:
        load_initial_data(fname, contig_data, checkm_data)

    first_comparison = pairwise_process(fnames[0], fnames[0], contig_data, checkm_data, overlap, 1)

    for idx, fname in enumerate(fnames[1:]):
        pairwise_process('current', fname, contig_data, checkm_data, overlap, idx+2)

    # Saving data
    df = checkm_data['current'].sort_values(['ID', 'Source'])
    df['NewID'] = range(1, len(df) + 1)

    outfile = os.path.join(outdir, 'binsources.tsv')
    df.to_csv(outfile, sep='\t', index=False)

    outfile = os.path.join(outdir, 'checkm.results.tsv')
    df = reconstruct_checkm_results(df, checkm_data)
    df.to_csv(outfile, sep='\t', index=False)

    save_contigs(df, outdir)

# -------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-checkmfiles',
                        required=True,
                        nargs='+',
                        help='Path to checkm result files list')

    parser.add_argument('-outdir',
                        required=True,
                        help='Path to outdir')

    parser.add_argument('-overlap',
                        required=False,
                        default=80,
                        type=int,
                        help='Minimum overlap (default 80)')

    args = parser.parse_args()

    checkmfiles = args.checkmfiles
    if len(checkmfiles) < 2:
        raise IOError("Please provide at least two checkM result files")

    process(checkmfiles, args.outdir, args.overlap)