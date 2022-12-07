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

    df_count["count_HTevents"] = df_count.apply(lambda x: count_HT(x,tree), axis = 1)
    df_count["spot_number"] = df_count.apply(lambda x: x["Family"].split("_")[0], axis = 1)
    df_count[["spot_number", "Family", "count_HTevents"]].to_csv(f"{outdir}/2-spot_pangenome/Final_spot_pangenome_HTg.tsv", sep = "\t", index = False)


def count_HT(row, tree):
    #considering one HTgene as a gain in a terminal branch >0.95 and excluding others terminals branchs with a last common ancestor with a gain >0.5
    #Also excluding HTgene gained at the last common ancestor of the phylogeny (meaning there are not part of the core genome just because of deletion event and not gain)
    list_ter_gain = []
    list_nonter_gain = []
    list_nonter = [x.name for x in tree.get_nonterminals() if x]
    list_ter = [x.name for x in tree.get_terminals() if x]

    for node, gain_value in row.items():
        if node in list_nonter and gain_value >= 0.5:
            list_nonter_gain.append(node)
        elif node in list_ter and gain_value >= 0.95:
            list_ter_gain.append(node)
    
    #Now if there are gain(s) of the gene we check how many differents event of horizontal transfer acquisition exists
    gain_event_count = 0
    list_already_counted = []
    d_gain = {}
    for nonter_node in list_nonter_gain:
        d_gain[nonter_node] = []
        for ter_node in list_ter_gain:
            if sorted(tree.find_clades(nonter_node))[0].is_parent_of(target = ter_node) and ter_node not in list_already_counted:
                list_already_counted.append(ter_node)
                d_gain[nonter_node].append(ter_node)

    gain_event_count = len(d_gain.keys())
    #adding all single terminal nodes which acquired the gene "by themselves"
    for ter_node in list_ter_gain:
        if ter_node not in list_already_counted:
            gain_event_count += 1
    
    return gain_event_count

    


def HTg_per_spot(outdir):
    #small function which create a tsv file with 5 columns which are as follow: spot_number, number_of_Htevents, number_of_accessory_genes_families, ref_left_c, ref_right_c
    #Also return two values, the number total of HTevents and the number total of spots
    df_htevents = pd.read_csv(f"{outdir}/2-spot_pangenome/Final_spot_pangenome_HTg.tsv", sep ="\t")
    df_ref_intervals = pd.read_csv(f"{outdir}/1-core_intervals/Ref_genome_intervals.tsv", sep = "\t") #useful to get the reference coordinates of the spot + the number total of spot
    n_total_spots = max(df_ref_intervals["interval_number"].tolist()) #Note first index is 0

    d_HTg_per_spot = {}
    n_total_HT_events= 0
    n_total_accessory_gene_families = 0
    for i in range (n_total_spots+1):
        if df_htevents[df_htevents["spot_number"] == str(i)].empty == False:
            d_HTg_per_spot[i] = {"spot_number": i,
                                "HTevents_number": sum(df_htevents[df_htevents["spot_number"] == str(i)]["count_HTevents"].tolist()),
                                "n_accessory_gene_families": len(df_htevents[df_htevents["spot_number"] == str(i)]),
                                "ref_genome_left_c": df_ref_intervals[df_ref_intervals["interval_number"] == i]["left_c"].values[0],
                                "ref_genome_right_c": df_ref_intervals[df_ref_intervals["interval_number"] == i]["right_c"].values[0]
                                }
            n_total_HT_events += sum(df_htevents[df_htevents["spot_number"] == str(i)]["count_HTevents"].tolist())
            n_total_accessory_gene_families += len(df_htevents[df_htevents["spot_number"] == str(i)])

        else :
            d_HTg_per_spot[i] = {"spot_number": i,
                                "HTevents_number": 0,
                                "n_accessory_gene_families": 0,
                                "ref_genome_left_c": df_ref_intervals[df_ref_intervals["interval_number"] == i]["left_c"].values[0],
                                "ref_genome_right_c": df_ref_intervals[df_ref_intervals["interval_number"] == i]["right_c"].values[0]
                                } 

    #adding a final row to store the total of HTevents in the whole analysis + total number of genes families
    d_HTg_per_spot[n_total_spots+1] = {"spot_number": "TOTAL HTevents/TOTAL accessory families",
                        "HTevents_number": n_total_HT_events ,
                        "n_accessory_gene_families": n_total_accessory_gene_families,
                        "ref_genome_left_c": "",
                        "ref_genome_right_c": ""
                        }
    
    
    pd.DataFrame.from_dict(d_HTg_per_spot, orient = "index").to_csv(f"{outdir}/3-HC_analysis/HTevents_per_spot.tsv", sep = "\t", index = False)

    return n_total_HT_events, n_total_spots,