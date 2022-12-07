import pandas as pd
import random
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# This script aims to determine which spot can be considered as a hotspot or a cold spot
# This following the protocol in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5635113/ publication, supplementary Fig 4



def HC_sim(n_HT_genes, m_spots, outdir, **kwargs):
    mode = kwargs.get("mode", "operon") # to determine which simulation we need to run
    n_runs = kwargs.get("n_runs", 1000)

    n = 0
    d_res_sim = {}
    for i in tqdm(range(n_runs), desc = "Running the simulation of HTgenes distribution"):
        d_res_sim[i] = {"max_HTgene_1spot": run_sim(n_HT_genes, m_spots, mode)}

    #now creating the histogram plot
    t95 = distribution_plot(d_res_sim, outdir)
    return t95


def run_sim(n_HT_genes, m_spots, mode):

    d_res = {} # to store in which spot are randomly associated X HT genes
    d_res = init_dict(d_res,m_spots)

    if mode =="operon":
        #Random localization of n_HT_genes*2/3
        for i in range(int(n_HT_genes*2/3/3)):
            d_res[choose_random_spot(m_spots)] += 3 #adding 3 genes to this spot
        
        for i in range(int(n_HT_genes*1/3)):
            d_res[choose_random_spot(m_spots)] += 1


    elif mode == "singleton":
        for i in range(int(n_HT_genes)):
            d_res[choose_random_spot(m_spots)] += 1

    #now getting the higher value of HT genes in a same interval and return it
    return max(d_res.values())

def init_dict(d, key_range):

    #small function which for i in range key_range create the key = i and value = 0 in the dict d

    for i in range(1, key_range+1):
        d[i] = 0
    
    return d

def choose_random_spot(max_number_spot):
    
    return random.randint(1,max_number_spot)

def distribution_plot(d_res_sim, outdir):
    #first need to create a dataframe
    df_sim = pd.DataFrame.from_dict(d_res_sim, orient = "index").reset_index()

    plt.hist(df_sim["max_HTgene_1spot"], color = "blue", edgecolor = "black", bins = int(180/5))
    plt.title("Distribution simulation of HTgenes random spot allocations")
    plt.xlabel("max_HTgene_1spot")
    plt.ylabel("n_simulation_run")
    plt.tight_layout()
    plt.savefig(f"{outdir}/3-HC_analysis/HC_simulations_distribution.png")

    with open(f"{outdir}/3-HC_analysis/HC_simulations_distribution.tsv", "w") as f:
        f.write(f"#T95% on this dataset = {np.percentile(df_sim['max_HTgene_1spot'].tolist(), 95)}\n")
    with open(f"{outdir}/3-HC_analysis/HC_simulations_distribution.tsv", "a") as f:
        df_sim.to_csv(f, sep = "\t", index = False)

    return np.percentile(df_sim['max_HTgene_1spot'].tolist(), 95)

def HC_spot(t95, outdir):

    #function which update 
    df_spot = pd.read_csv(f"{outdir}/3-HC_analysis/HTevents_per_spot.tsv", sep = "\t") 
    pd.to_numeric(df_spot["HTevents_number"])
    df_spot["HC_spot"] = df_spot.apply(lambda x: "Hotspot" if x["HTevents_number"] >= t95 else ("Coldspot" if x["n_accessory_gene_families"] > 0 else "Empty_spot"), axis = 1)
    index_sum = df_spot[df_spot["spot_number"] == "TOTAL HTevents/TOTAL accessory families"].index[0]
    df_spot.loc[index_sum,"HC_spot"] = ""
    df_spot.to_csv(f"{outdir}/3-HC_analysis/HC_spots_pangenome.tsv", sep = "\t", index = False)

    #Also produce a summary file
    with open(f"{outdir}/3-HC_analysis/HC_summary.txt", "w") as f:
        f.write(f"N_accessory_genes_whole_pangenome:\t{df_spot.loc[index_sum,'n_accessory_gene_families']}\n")
        f.write(f"N_HTevents_whole_pangenome:\t{df_spot.loc[index_sum,'HTevents_number']}\n")
        f.write(f"T95_percentile:\t{t95}\n")
        f.write(f"N_total_spots:\t{index_sum}\n")
        f.write(f"N_Hotspot(s):\t{len(df_spot[df_spot['HC_spot'] == 'Hotspot'])} ({len(df_spot[df_spot['HC_spot'] == 'Hotspot'])/index_sum*100}% all spots or {len(df_spot[df_spot['HC_spot'] == 'Hotspot'])/(len(df_spot[df_spot['HC_spot'] == 'Hotspot'])+len(df_spot[df_spot['HC_spot'] == 'Coldspot']))*100}% non empty spots)\n")
        f.write(f"N_Coldspot(s):\t{len(df_spot[df_spot['HC_spot'] == 'Coldspot'])} ({len(df_spot[df_spot['HC_spot'] == 'Coldspot'])/index_sum*100}% all spots or {len(df_spot[df_spot['HC_spot'] == 'Coldspot'])/(len(df_spot[df_spot['HC_spot'] == 'Hotspot'])+len(df_spot[df_spot['HC_spot'] == 'Coldspot']))*100}% non empty spots)\n")
        f.write(f"N_Empty_spot(s):\t{len(df_spot[df_spot['HC_spot'] == 'Empty_spot'])} ({len(df_spot[df_spot['HC_spot'] == 'Empty_spot'])/index_sum*100}%)\n")
        f.write(f"HTg50:\t{HTg50(df_spot)}\n")


def HTg50(df):
    #return the minimal number of spot to contains 50% of the HTgenes
    total_HTevents = df.loc[df[df["spot_number"] == "TOTAL HTevents/TOTAL accessory families"].index[0], "HTevents_number"]
    df = df[df["spot_number"] != "TOTAL HTevents/TOTAL accessory families"].sort_values(by= "HTevents_number", ascending = False).reset_index(drop =True)

    count = 0
    count_htevent = 0
    half_htevent = int(total_HTevents)/2

    while count_htevent < half_htevent:
        count_htevent += df.loc[count,"HTevents_number"]
        print(count_htevent, "/", half_htevent)
        count += 1
        print("Count",count)
    
    return count





#just for testing
if __name__ == "__main__":
    HC_sim(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

    #python script/HCspot.py 200 1722 results/2022-12-02_TATdyn/3-HC_analysis