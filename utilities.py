import os
import glob
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from Bio import Phylo
import multiprocessing
from multiprocessing import  Pool


def check_dir(outdir, date, **kwargs):
    TAT_defense = kwargs.get("TAT_defense", False)

    list_f = ["1-core_intervals", "2-spot_pangenome", "3-HC_analysis", "4-Spot_dissimilarity", "tmp", "tmp/mmseqs"]
    if TAT_defense == True:
        list_f.append("5-optional_defense-finder")

    list_folders2create = []
    list_folders2create.append(outdir)
    list_folders2create.append(f"{outdir}/{date}_TATdyn")
    for i in list_f:
        list_folders2create.append(f"{outdir}/{date}_TATdyn/{i}")

    for i in list_folders2create:
        try :
            os.mkdir(i)
        except FileExistsError:
            continue


def ref_genome(input_folder, user_ref):

    # reading the namefiles of each genome in the input folder
    list_file = glob.glob(f"{input_folder}/gff3/*.gff")
    list_genome = []
    for i in list_file:
        list_genome.append(i.rsplit("/",1)[-1].rsplit(".",1)[0])
    
    if user_ref != None :
        if user_ref in list_genome:
            ref = user_ref
        else :
            raise Exception(f"{user_ref} is not part of the genomes given")
    
    else :
        ref = list_genome[0]
    
    #checking if the given genome genome is fully assembled
    df_ref = pd.read_csv(f"{input_folder}/gff3/{ref}.gff", sep = "\t", comment = "#", names = ["contig", "source", "type", "left_c", "right_c", ".", "strand", "0", "id"])
    n_contigs = len(df_ref["contig"].drop_duplicates().tolist())
    if n_contigs == 1:
        return ref, list_genome

    else :
        raise Exception(f"{ref} is not a fully assembled genome, please provide a fully assembled genome for the reference (detected contig : {n_contigs} contigs") 
    

def combine_input(input_folder, outdir, list_genomes):

    #return two var: a multi index SeqIO object
    #and a concatanated dataframe with all dataframe formated from the user input (note: it will add a column with the genome name to facilitate the use later)

    #first creating the multi gff dataframe
    list_df2concatenate = []
    
    for i in list_genomes :
        df_tmp = pd.read_csv(f"{input_folder}/gff3/{i}.gff", sep = "\t", comment = "#", names = ["contig", "source", "type", "left_c", "right_c", ".", "strand", "0", "id"])
        df_tmp["genome"] = i
        list_df2concatenate.append(df_tmp)
    
    df_all_gff = pd.concat(list_df2concatenate)
    df_all_gff["id"] = df_all_gff["id"].apply(lambda x: x.split(";",1)[0].split("=")[-1])
    df_all_gff.to_csv(f"{outdir}/tmp/concat_gff.gff", sep = "\t", index = False)

    # now creating a multi index SeqIO object for 
    str_cat_prt = "cat "
    str_cat_fna = "cat "

    for genome in list_genomes:
        str_cat_fna += f"{input_folder}/Replicons/{genome}.fna "
        str_cat_prt += f"{input_folder}/Proteins/{genome}.prt "
    
    str_cat_prt += f"> {outdir}/tmp/cat_proteins.prt"
    str_cat_fna += f"> {outdir}/tmp/cat_replicons.fna"
    os.system(str_cat_prt)
    os.system(str_cat_fna)


def threads2use(n):
    executor = ThreadPoolExecutor()

    if int(n) != 0 and int(n) <= executor._max_workers:
        return n
    else:
        return executor._max_workers


def core2use(n):
    if int(n) != 0 and int(n) <= multiprocessing.cpu_count():
        return n
    else :
        return multiprocessing.cpu_count()

def name_tree(treefile, outdir):
    tree = Phylo.read(treefile, "newick")
    n = 1
    for node in tree.get_nonterminals():
        node.name = f"node{n}"
        n += 1
    
    Phylo.write(tree, f"{outdir}/2-spot_pangenome/tree_renamed.nwk", 'newick')


def multicore_pd_apply(df, func, **kwargs):
    #small function to parralelize pandas apply functions
    n_core = kwargs.get("n_core", 1)

    df_split = np.array_split(df,n_core)
    pool = Pool(n_core)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


def combine_results(outdir, **kwargs):
    user_genes = kwargs.get("user_gene", None)
    TAT_defense = kwargs.get("TAT_defense", False)
    #small functions which will combine the results from Hot/Coldspot analysis with the dissimilarity index for each spot
    # Also will add a column with the localization of each user genes if given
    #loading both results dataframe

    df_HC = pd.read_csv(f"{outdir}/3-HC_analysis/HC_spots_pangenome.tsv", sep = "\t")
    df_dissimilarity = pd.read_csv(f"{outdir}/4-Spot_dissimilarity/Final_spot_dissimilarity_index.tsv", sep = "\t")

    #reordering the dataframe + renaming the columns
    df_HC = df_HC.rename(columns = {"spot_number": "spot", "HTevents_number": "HTg", "n_accessory_gene_families": "n_acc", "ref_genome_left_c": "ref_left_c", "ref_genome_right_c": "ref_right_c"})
    df_HC = df_HC.loc[:,["spot", "ref_left_c", "ref_right_c", "n_acc", "HTg", "HC_spot"]].iloc[:-1]
    df_HC["spot"] = df_HC["spot"].astype("int")
    df_dissimilarity["spot"] = df_dissimilarity["spot"].astype("int")
    df_res = pd.merge(df_HC, df_dissimilarity, on = "spot", how = "left")
    df_res[["ref_left_c", "ref_right_c"]] = df_res[["ref_left_c", "ref_right_c"]].astype("int")

    #in case the user provide a file with the coordinates of elements to localize, add a new column with the id of the element within the corresponding spot
    if user_genes != None:
        df_user = pd.read_csv(user_genes, sep = "\t", comment = "#", names = ["id", "left_c", "right_c"])
        df_user[["left_c", "right_c"]] = df_user[["left_c", "right_c"]].astype("int")
        df_res["user_id"] = [[] for _ in range(len(df_res))]
        for elements_row in df_user.iterrows():
            if df_res[(df_res["ref_left_c"] <= int(elements_row[1]["left_c"])) & (df_res["ref_right_c"] >= int(elements_row[1]["right_c"]))].empty == False:
                tmp_idx = df_res[(df_res["ref_left_c"] <= int(elements_row[1]["left_c"])) & (df_res["ref_right_c"] >= int(elements_row[1]["right_c"]))].index[0]
                df_res.loc[tmp_idx, "user_id"].append(elements_row[1]['id'])

            else : 
                print(f"User input element {elements_row[1]['id']}, could not be placed within the spot pangenome generated (part of core genome ?)")

        df_res["user_id"] = df_res.apply(lambda x: (",").join(x["user_id"]), axis = 1)

    if TAT_defense == True:
        df_defense = pd.read_csv(f"{outdir}/5-optional_defense-finder/1bis_spot_pangenome_defense_system_numbers.tsv", sep= "\t")
        for defense_spot in df_defense.iterrows():
            str2add = ""
            total_defense = 0
            for defense_type in df_defense.loc[:,df_defense.columns != "Spot_number"].columns.tolist():
                total_defense += defense_spot[1][defense_type]
                if defense_spot[1][defense_type] != 0:
                    str2add += f"{defense_spot[1][defense_type]} {defense_type} ~ "

            df_res.loc[defense_spot[1]["Spot_number"],"Defense_systems"] = f"{total_defense} defense systems: {str2add}"

    df_res.to_csv(f"{outdir}/all_results_combined.tsv", sep = "\t", index = False)
    if user_genes != None :
        df_res[df_res["user_id"] != ""].to_csv(f"{outdir}/all_results_with_user_id_spots_only.tsv", sep = "\t", index = False)
    df_filtered = df_res[df_res["HC_spot"] != "Empty_spot"]
    df_filtered = df_filtered[(df_filtered["βsor"] != 0) & (df_filtered["βsor"].notnull()) & (df_filtered["βsim"] != "div0")].sort_values(by ="HC_spot")
    df_filtered.to_csv(f"{outdir}/all_results_combined_filtered.tsv", sep = "\t", index = False)
    
