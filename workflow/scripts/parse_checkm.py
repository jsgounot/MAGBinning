# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2021-07-02 15:50:09
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-07-21 15:35:53

import os, glob, re
import pandas as pd

bname = os.path.basename
dname = os.path.dirname

def touch(fname):
    with open(fname, "w") as f:
        pass

def load_eval(evalFile) :
    data = []
    
    # Wanted column index
    columns = {
        "ID" : 1,
        "Completeness" : 6,
        "Contamination" : 7,
        "Size" : 9,
        "#Contigs" : 11,
        "Contig_N50" : 14,
        "Longest_contig" : 18
        }
    
    with open(evalFile) as f :
        for line in f :
            line = re.split(" \s+", line)
            if len(line) == 1 or line[1] == "Bin Id" : continue
            data.append({column : line[idx] for column, idx in columns.items()})
                                       
    assert len(evalFile)
                                                
    df = pd.DataFrame(data)
    
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

fname = snakemake.input[0]
outfile = snakemake.output[0]

if os.path.basename(fname) == "empty_file.tmp":
    touch(outfile)
    exit(0)
else:
    df = load_eval(fname)
    df.to_csv(outfile, sep="\t")