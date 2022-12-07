import os
import glob
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from Bio import Phylo

def check_dir(outdir, date):

    list_f = ["1-core_intervals", "2-spot_pangenome", "tmp", "tmp/mmseqs"]

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


def name_tree(treefile, outdir):
    tree = Phylo.read(treefile, "newick")
    n = 1
    for node in tree.get_nonterminals():
        node.name = f"node{n}"
        n += 1
    
    Phylo.write(tree, f"{outdir}/2-spot_pangenome/tree_renamed.nwk", 'newick')