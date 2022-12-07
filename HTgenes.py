import pandas as pd
from Bio import Phylo
import os
from datetime import datetime

def run_Count(outdir, **kwargs):
    count_path = kwargs.get("count_path", os.path.expanduser("~/count/Count/Count.jar"))

    #First creating the birth/death/duplication model using ML 
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Currently running Count to generate the birth/death/duplication model")
    os.system(f"java -cp {count_path} -Xmx2048M ca.umontreal.iro.evolution.genecontent.ML {outdir}/2-spot_pangenome/tree_renamed.nwk {outdir}/2-spot_pangenome/whole_pangenome_summary.table > {outdir}/2-spot_pangenome/whole_pangenome_bdd.rates")
    print(f"[{datetime.now().strftime('%H:%M:%S')}] Currently running Count with the pangenome to detect HTgs")
    os.system(f"java -cp {count_path} -Xmx2048M ca.umontreal.iro.evolution.genecontent.Posteriors {outdir}/2-spot_pangenome/tree_renamed.nwk {outdir}/2-spot_pangenome/spot_pangenome_summary_formated.table {outdir}/2-spot_pangenome/whole_pangenome_bdd.rates > {outdir}/2-spot_pangenome/count_spot_pangenome.xml")


def determine_HTgenes(outdir):
    #function which determine which gene from a given pangenome is considered as Horizontal Transfer genes (HTg)
    # This using count and the birth - loss - duplication model

    df_count = pd.read_csv(f"{outdir}/2-spot_pangenome/count_spot_pangenome.xml", sep = "\t", comment = "#")
    #loading the tree so we can't determine if the last common ancestor is the one which gain the gene or not
    tree = Phylo.read(f"{outdir}/2-spot_pangenome/tree_renamed.nwk", "newick")

    list_columns2keep = [x for x in df_count.columns.tolist() if ":gain" in x]
    list_columns2keep.insert(0,"Family")
    df_count = df_count[list_columns2keep]
    list_columns2rename = [x.replace(":gain","") for x in list_columns2keep]

    d_rename = {}
    n = 0
    for i in list_columns2keep:
        d_rename[i] = list_columns2rename[n]
        n += 1
    df_count = df_count.rename(columns = d_rename)
    list_nonter = [x.name for x in tree.get_nonterminals()]
    print(list_nonter)
    print(tree.root.name)
    #df_count["count_HTg"] = df_count.apply(lambda x: count_HT(x,tree), axis = 1)


def count_HT(row, tree):
    #considering one HTgene as a gain in a terminal branch >0.95 and excluding others terminals branchs with a last common ancestor with a gain >0.5
    #Also excluding HTgene gained at the last common ancestor of the phylogeny (meaning there are not part of the core genome just because of deletion event and not gain)
    n_HT_events = 0
    list_nonter = [x.name for x in tree.get_nonterminals()]
    list_nonter = list_nonter.remove(tree.root.name)
    list_ter = [x.name for x in tree.get_terminals()]
    print(list_nonter)



