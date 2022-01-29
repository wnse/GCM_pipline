import os
import argparse
import scipy
import gower
import pandas as pd
import numpy as np


def make_cgMLST_tree(alleles_tsv):
    tree_file = gower.gower_matrix(df.astype(np.float32))
    return tree_file

def merge_alleles_tsv(db_alleles, chewbbaca_alleles, sampleName):
    df_ref = pd.read_csv(db_alleles, sep='\t', index_col=0)
    df_query = pd.read_csv(chewbbaca_alleles, sep='\t', index_col=0)
    df_query.index = df_query.index.str.replace(r'\.fasta$', '', regex=True)
    df_query = df_query.apply(lambda a: a.where(a.astype(str).str.isnumeric()).fillna(0), axis=1)
    core_gene = np.intersect1d(df_query.columns, df_ref.columns)
    df_merge = pd.concat([df_ref[core_gene], df_query[core_gene]])

    return df_merge

