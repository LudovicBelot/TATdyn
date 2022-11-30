import os
import glob
import pandas as pd


def check_dir(outdir, date):

    list_f = ["1-core_intervals", "tmp"]

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
        return ref

    else :
        raise Exception(f"{ref} is not a fully assembled genome, please provide a fully assembled genome for the reference (detected contig : {n_contigs} contigs") 