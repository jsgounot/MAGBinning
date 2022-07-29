# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-07-25 16:00:51
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-07-29 10:30:50

'''
Checkm2 does not return the ['Size', '#Contigs', 'Contig_N50', 'Longest_contig'] like checkm.
So we need to compute those ourself here
'''

import os, glob
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

def touch(fname):
    with open(fname, "w") as f:
        pass

def getn50(l):
    l = sorted(l)[::-1]
    s = sum(l) / 2
    c = 0
    for idx, value in enumerate(l):
        c += value
        if c >= s:
            return value

def fasta_info(fname):
    rec = []
    l = []

    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line[0] == '>':
                if rec: l.append(len(''.join(rec)))
                rec = []
            else:
                rec.append(line)

    if rec: l.append(len(''.join(rec)))

    return {
        'ID': os.path.basename(fname),
        'Size' : sum(l),
        '#Contigs': len(l),
        'Contig_N50': getn50(l),
        'Longest_contig': max(l)
    }
def retrieve_finfo(binlist):
    with open(binlist) as f:
        return pd.DataFrame(fasta_info(fname.strip())
        for fname in f)

def get_fa_id(binlist):
    with open(binlist) as f:
        d = {}
        for fname in f :
            fname = fname.strip()
            if not fname: continue
            bname = os.path.basename(fname)
            d[os.path.splitext(bname)[0]] = bname
        return d

def load_quality_report(quality_report) :

    df = pd.read_csv(quality_report, sep='\t')

    df.columns = [{'Name': 'ID'}.get(column, column)
        for column in df.columns]

    df = df[['ID', 'Completeness', 'Contamination']]

    def set_completeness(row) :
        completeness = row["Completeness"]
        contamination = row["Contamination"]
        if completeness > 90 and contamination < 5 :
            return "HIGH"
        elif completeness >= 50 and contamination <= 10 :
            return "MEDIUM"
        else :
            return "LOW"
                                                                                                
    df["Contamination"] = df["Contamination"].astype(float)
    df["Completeness"] = df["Completeness"].astype(float)                                                                                            
    df["CheckMStatus"] = df.apply(lambda row : set_completeness(row), axis=1)

    return df

report = snakemake.input['report']
binlist = snakemake.input['binlist']
outfile = snakemake.output[0]

if os.path.basename(report) == "empty_file.tmp":
    touch(outfile)
    exit(0)
else:
    df = load_quality_report(report)
    gfi = get_fa_id(binlist)
    df['ID'] = df['ID'].apply(lambda mid: gfi[mid])
    sdf = retrieve_finfo(binlist)
    df = df.merge(sdf, on='ID', how='left')

    if len(sdf) != len(df):
        raise Exception(f'Dataframe are not the same size ! {report} {binlist}')

    df.to_csv(outfile, sep="\t")