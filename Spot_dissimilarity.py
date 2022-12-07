import pandas as pd
from tqdm import tqdm





def prepare_data(outdir):
    #function which uses the previous results to create a tsv files with all informations needed for spot turneover/nestedness calcul

    df_table_spot_mmseqs = pd.read_csv(f"{outdir}/2-spot_pangenome/spot_pangenome_summary.table", sep = "\t").rename(columns = {"Unnamed: 0":"id"})
    df_table_spot_mmseqs["spot"] = df_table_spot_mmseqs.apply(lambda x: x["id"].split("_")[0], axis = 1)

    #creating an index genome name {g1:"genome_name1", ...} to facilitate the following script writing
    index_genome = 0
    d_genome_index = {}
    d_genome_index4file = {}
    for i in df_table_spot_mmseqs.columns.tolist():
        if i != "id" and i != "spot":
            index_genome+=1
            d_genome_index4file[index_genome] = {"Real_genome_name": i, "idx" : f"g{index_genome}"}
            d_genome_index[f"g{index_genome}"] = i

    pd.DataFrame.from_dict(d_genome_index4file, orient = "index").to_csv(f"{outdir}/4-Spot_dissimilarity/genome_index.idx", sep ="\t", index = False)
    
    list_genome_combinaisons = generate_combinaison(index_genome)
    
    d_res = {}
    test = 0
    for spot in tqdm(df_table_spot_mmseqs["spot"].tolist(), desc = "Preparing each spot data for dissimilarity calcul"):
        d_res[spot] = {}
        for g, genome_name in d_genome_index.items():
            d_res[spot][f"{g}_total_accessory"] = count_accessory_genes(genome_name, spot, df_table_spot_mmseqs)
        for genome_combinaison in list_genome_combinaisons:
            d_res[spot][f"g{genome_combinaison[0]}_g{genome_combinaison[1]}_common_accessory"] = count_common_accessory_genes(d_genome_index[f"g{genome_combinaison[0]}"], d_genome_index[f"g{genome_combinaison[1]}"], spot, df_table_spot_mmseqs)

    df_res = pd.DataFrame.from_dict(d_res, orient = "index")
    df_res["ST"] = df_res.apply(lambda x: total_accessory(x), axis = 1)

    df_res.to_csv(f"{outdir}/4-Spot_dissimilarity/data_prepared4calcul.tsv", sep ="\t", index = False)



def generate_combinaison(j):
    #generate combinaisons of genomes for which , considering a number of j genomes , i < j (example for j = 3, it exists 3 combinaisons => g1,g2    g1,g3   g2,g3)
    #return a list of tuples
    list_of_combinaisons = []

    last_i = 1
    while last_i<j:
        for i in range(1,j+1):
            if i > last_i:
                list_of_combinaisons.append((last_i,i))
        last_i += 1
    
    return list_of_combinaisons


def count_accessory_genes(genome_name, spot_number, df_mmseqs):
    # function which for a given genome g count the number of accessory genes at the given spot
    tmp_list = df_mmseqs[df_mmseqs["spot"] == spot_number][genome_name].tolist()
    if "?" in tmp_list : #couldn't recreate the interval in this genome
        return "?"
    else :
        tmp_list = [int(x) for x in tmp_list if x]
        return sum(tmp_list)


def count_common_accessory_genes(genome1, genome2, spot_number, df_mmseqs):
    # function which for two given genomes g1,g2 count the number of accessory genes in common between both genomes at the given spot
    df = df_mmseqs[df_mmseqs["spot"] == spot_number][[genome1,genome2]]

    #checking if we could recreate the interval in both genome first
    if "?" in df[genome1].tolist() or "?" in df[genome2].tolist():
        return "?"
    else :
        count = 0
        for gene_family in df.iterrows():
            if int(gene_family[1][genome1]) > 0 and int(gene_family[1][genome2]) > 0:
                count += min(int(gene_family[1][genome1]), int(gene_family[1][genome2]))
        
        return count


def total_accessory(row):
    count = 0
    for i in row.values:
        if i != "?":
            count += int(i)
    
    return count