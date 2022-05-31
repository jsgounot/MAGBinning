# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 15:36:26
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-05-31 14:57:10

import os
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

tRNAStrict = set([
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        ])

"""
tRNATCount : Total number of tRNA including undefined results
tRNACount : Number of tRNA without undefined results but with redundancy
tRNAUnique : Unique number of tRNA including non-strict tRNA
tRNAUniqueStrict : Unique number of tRNA considering only strict ones (the one which should be considered for MIMAG analysis)
"""

def touch(fname):
    with open(fname, "w") as f:
        pass

# ------------------------------------------------------

inputs = list(snakemake.input)
outfile = snakemake.output[0]
if os.path.basename(inputs[0]) == 'binlist.txt' :
    # Case where no bins where found
    touch(outfile)

else :
    trnadata = []

    for fname in sorted(inputs) :

        kind = os.path.basename(os.path.dirname(fname))
        sidx = os.path.basename(fname).split('.')[1]
        mID  = os.path.basename(snakemake.params[0][sidx])

        df = pd.read_csv(fname, sep="\t", usecols=[0,4], names=["name", "RNAType"], skiprows=3)
        totalu = len(df)

        df = df[df["RNAType"] != "Undet"] # pseudo rna
        df = df[~ df['RNAType'].isnull()]

        if df.empty :
            total = nunique = nunique_strict = 0

        else :
            total, nunique = len(df), df["RNAType"].nunique()
            nunique_strict = len(set(df["RNAType"].apply(lambda name : name.upper())) & tRNAStrict)
        
        trnadata.append(
            {"ID" : mID,
            "tRNAKind" : kind,
            "tRNACount" : total, 
            "tRNAUnique" : nunique, 
            "tRNAUniqueStrict" : nunique_strict, 
            "tRNATCount" : totalu}
            )

    df = pd.DataFrame(trnadata)
    df.to_csv(outfile, sep="\t")