# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-10-21 22:14:19
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-31 10:54:46

"""
This is a new version of the metaWrap script `binning_refiner.py` (1), itself coming from the 
binning_refiner script (2). Note that the original author published an update of the metaWrap 
script (a R script) which include some of the modification I add here.

What this script doing is making new bins based on previously done bins. In other word,
if you have binner A and binner B, with bins A1, A2, A3 and B1, B2. This script will look at 
all possible combinations (A1, B1), (A1, B2), etc. For each combinations, we look at all contigs
found both in A1 and B1, and if the sum of all contigs length is higher than minsize, a new bins
with contigs found in both of them (THE INTERSECTION) will be formed based on these contigs.

NOTE THAT BINNING MUST COME FROM THE SAME ASSEMBLY SINCE MERGING IS BASED ON CONTIGS NAME

Comparison with the metaWrap script:
- Clean the messy code, use pandas instead of a crazy amount of dictionaries
- Do not write temporary fasta along the way
- Allow unlimited number of bins
- Allow multiple cores
- Allow to decompose all combinations in one call
- Check that fasta sequence are the same based on sequences hash
- Should be in general much faster

Input variation:
- Instead of a directory, this script takes a list of files containing path of bins
  which can be easily created with a simple `ls -d *.fasta > binlist.txt'

I keep the weird size threshold use in the original script to have the same result compared
to previous script.

(1) https://github.com/bxlab/metaWRAP/blob/master/bin/metawrap-scripts/binning_refiner.py
(2) https://github.com/songweizhi/Binning_refiner
"""

import argparse, json, os, string
import pandas as pd
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
from Bio import SeqIO

def extract_fnames(binlistfile):
    # Open a file containing path
    # Return paths as list and check they exist

    fnames = []
    with open(binlistfile) as f:
        for line in f:
            line = line.strip()
            if not os.path.isfile(line): 
                raise OSError('File not found: ' + line)
            fnames.append(line)
    return fnames

def read_fasta(fname):
    with open(fname) as f:
        for record in SeqIO.parse(f, 'fasta'):
            yield fname, record.id, len(record.seq), hash(str(record.seq))

def process_fasta(fname):
    # Return a dataframe with each row being a contig
    # with the file path where it comes from, the contig ID,
    # the sequence length and the sequence hash value

    return pd.DataFrame([
        {
            'filename': fname,
            'seqid': seqid,
            'seqlen': seqlen,
            'seqhash': seqhash
        }
        for fname, seqid, seqlen, seqhash in read_fasta(fname)
    ])

def process_fasta_ncore(fnames, outfile, ncore=1):
    # Multiprocess a list of fasta file paths
    # For each file, run process_fasta and return
    # a dataframe containing info from all files

    pool = Pool(ncore)
    args = [(fname,) for fname in fnames]

    try :
        data = pool.starmap(process_fasta, args)
        pool.close()
        pool.join()

    except KeyboardInterrupt :
        print("Interrupt childs process")
        pool.terminate()
        raise KeyboardInterrupt()

    return pd.concat(data)

def extract_records(fname):
    with open(fname) as f:
        for record in SeqIO.parse(f, 'fasta'):
            yield record

def process_group(df, fnacols, outdir, minsize):
    # We create new bins for each group containing contigs found in the same bins in all binners
    
    groupidx = 0
    fnames = []

    for group, sdf in df.groupby(fnacols):
        # We check the size
        if sdf['seqlen'].sum() < minsize :
            continue

        groupidx += 1
        binfname = sdf[fnacols[0]].iloc[0]
        contigs  = set(sdf['contig'])

        outfile = os.path.join(outdir, 'refined_bin_' + str(groupidx) + '.fa')
        records = [record for record in extract_records(binfname) if record.id in contigs]
        assert len(records) == len(contigs)
        SeqIO.write(records, outfile, 'fasta')
        fnames.append(str(Path(outfile).absolute()))

    print (str(groupidx) + ' refined bins processed in ' + outdir)

    outfile = os.path.join(outdir, 'binsinfo.tsv')
    df.to_csv(outfile, sep='\t')

    outfile = os.path.join(outdir, 'binlist.tsv')
    with open(outfile, 'w') as f:
        f.write('\n'.join(fnames) + '\n')

def process(binlists, outdir, minsize=524288, decompose=True, ncore=1):
    # Run the main process 

    # 1. Read the fasta files contained in each binlist and generate a dataframe
    # merged on contig ID. The idea is too see for each contig, if we can find it in another bins
    # the file (the bin) where it comes from
    # the sequence length and its hash value

    df = None
    for idx, binlist in enumerate(binlists):
        fastas = extract_fnames(binlist)
        sdf = process_fasta_ncore(fastas, ncore)

        # Order and rename columns
        sdf = sdf[['seqid', 'filename', 'seqlen', 'seqhash']]
        sdf.columns = ['contig', 'filename_' + str(idx), 'seqlen_' + str(idx), 'seqhash_' + str(idx)]

        # We merge contigs found in all binners
        if df is None: df = sdf
        else: df = df.merge(sdf, on='contig', how='outer')

    # define the len. We use here min function to define a non NaN number
    lencols = ['seqlen_' + str(idx) for idx in range(len(binlists))]
    df['seqlen'] = df[lencols].min(axis=1)

    # We check that contigs are the same using the sequence hash
    hashcols = ['seqhash_' + str(idx) for idx in range(len(binlists))]
    df['hunique'] = df[hashcols].nunique(axis=1) == 1
    sdf = df[~ df['hunique']]
    if not sdf.empty:
        raise Exception('Some contig sequence from different binners does not look the same !! \n ' + sdf)

    if decompose:
        # We give a letter instead of a number for clarity
        charidx = {idx: char for idx, char in enumerate(string.ascii_uppercase)}

        # We generate group for each combinations
        used = ''.join(str(idx) for idx in range(len(binlists)))
        groups = (group for gsize in range(2, len(used) + 1) for group in combinations(used, gsize))
        dnames = {}

        for group in groups:           
            name = ''.join(charidx[int(idx)] for idx in group)
            suboutdir = os.path.join(outdir, name)
            os.makedirs(suboutdir, exist_ok=True)

            fnacols = ['filename_' + idx for idx in group]
            process_group(df, fnacols, suboutdir, minsize)

            dnames[name] = [binlists[int(idx)] for idx in group]

        outfile = os.path.join(outdir, 'groupname.json')
        with open(outfile, 'w') as f:
            json.dump(dnames, f, indent=4)

    else:
        # We run the refinment 
        fnacols = ['filename_' + str(idx) for idx in range(len(binlists))]
        process_group(df, fnacols, outdir, minsize)

# -------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-binlists',
                        required=True,
                        nargs='+',
                        help='Path to bin files list')

    parser.add_argument('-outdir',
                        required=True,
                        help='Path to outdir')

    parser.add_argument('-minsize',
                        required=False,
                        default=524288,
                        type=int,
                        help='Minimum size for a combination to be consider (default 524288)')

    parser.add_argument('-decompose',
                        action='store_true',
                        help='Not only proceed all binners at once, but decompose all subbinner combinations')

    parser.add_argument('-threads',
                        required=False,
                        default=1,
                        type=int,
                        help='Numbe of threads to use (default 1)')

    args = parser.parse_args()

    binlists = args.binlists
    if len(binlists) < 2:
        raise IOError("Please provide at least two binlist file (one for each binner)")
    
    outdir = args.outdir
    if not os.path.isdir(outdir):
        raise IOError('Outdir not found: ' + outdir)

    process(binlists, outdir, args.minsize, args.decompose, args.threads)
    