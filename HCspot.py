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
    distribution_plot(d_res_sim, outdir)


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
    plt.savefig(f"{outdir}/HC_simulations_distribution.png")

    with open(f"{outdir}/HC_simulations_distribution.tsv", "w") as f:
        f.write(f"#T95% on this dataset = {np.percentile(df_sim['max_HTgene_1spot'].tolist(), 95)}\n")
    with open(f"{outdir}/HC_simulations_distribution.tsv", "a") as f:
        df_sim.to_csv(f, sep = "\t", index = False)

def HC_pangenome(outdir):

    #function which determine which spot can be considered as an Hotspot or not (based on the simulation)
    df_spot = pd.read_csv(f"{outdir}/2-spot_pangenome/HTgenes_spot.tsv", sep = "\t")

    with open(f"{outdir}/3-HC_analysis/HC_simulations_distribution.tsv", "r") as f:
        for line in f:
            T95_percentile = int(line.split("#T95% on this dataset = ")[1].strip())
            break
    
    pd.to_numeric(df_spot["n_HTgenes"])
    df_spot["HC_spot"] = df_spot.apply(lambda x: "Hotspot" if x["n_HTgenes"] >= T95_percentile else "Coldspot")
    df_spot.to_csv(f"{outdir}/2-spot_pangenome/HTgenes_spot.tsv", sep = "\t", index = False)




#just for testing
if __name__ == "__main__":
    HC_sim(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])

    #python script/HCspot.py 200 1722 results/2022-12-02_TATdyn/3-HC_analysis