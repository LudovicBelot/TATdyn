import pandas as pd
from tqdm import tqdm






def ref_interval(ref_genome, core_file, outdir):
    
    #First we divide the reference genome in intervals based on their core genome 
    df_core = pd.read_csv(core_file, sep = "\t", names = ["gene_name", "core_family", "genome_name", "contig", "gene_number", "left_coordinate", "right_coordinate", "strand"], skiprows = 1)   
    df_core_ref = df_core[df_core["genome_name"] == ref_genome].sort_values(by = "left_coordinate").reset_index(drop= True)

    d_intervals = {}
    for row_index in df_core_ref.index.tolist():
        if row_index + 1 < len(df_core_ref):
            d_intervals[row_index] = {  "interval_number": row_index,
                                        "left_core": df_core_ref.loc[row_index, "gene_name"],
                                        "left_core_family": df_core_ref.loc[row_index, "core_family"],
                                        "right_core": df_core_ref.loc[row_index+1, "gene_name"],
                                        "right_core_family": df_core_ref.loc[row_index+1, "core_family"],
                                        "left_c": df_core_ref.loc[row_index, "right_coordinate"] ,
                                        "right_c": df_core_ref.loc[row_index+1, "left_coordinate"]
                                        }

        elif row_index + 1 == len(df_core_ref):
            d_intervals[row_index] = {  "interval_number": row_index,
                                        "left_core": df_core_ref.loc[row_index, "gene_name"],
                                        "left_core_family": df_core_ref.loc[row_index, "core_family"],
                                        "right_core": df_core_ref.loc[0, "gene_name"],
                                        "right_core_family": df_core_ref.loc[0, "core_family"],
                                        "left_c": df_core_ref.loc[row_index, "right_coordinate"] ,
                                        "right_c": df_core_ref.loc[0, "left_coordinate"]
                                        }
    
    df_intervals = pd.DataFrame.from_dict(d_intervals, orient = "index")

    df_intervals.to_csv(f"{outdir}/1-core_intervals/Ref_{ref_genome}_intervals.tsv", sep = "\t", index = False)




def check_interval(ref_genome, list_genomes, core_file, outdir):
    # Now based on the intervals generated with the reference genome, we try to rebuild the same interval in the others genomes
    # For that we need to check two things to accept the interval
    # First: in case the genome isn't fully assembled, both core genes used to delimit the interval need to be located on the same contig
    # Second : To make sure there are no important genomic rearrangement, for that we check if there are no others core gene within the same interval in this genome

    #it will return a dataframe with row = intervals , columns = genomes ====> each cell with either : accepted (we're able to rebuild the interval in this genome),
    # rejected (differents contigs) or rejected (rearrangments)
    # last columns count the number of accepted genomes in which we can rebuild the interval

    df_ref = pd.read_csv(f"{outdir}/1-core_intervals/Ref_{ref_genome}_intervals.tsv", sep = "\t")
    df_core = pd.read_csv(core_file, sep = "\t", names = ["gene_name", "core_family", "genome_name", "contig", "gene_number", "left_coordinate", "right_coordinate", "strand"], skiprows = 1)   
    df_gff = pd.read_csv(f"{outdir}/tmp/concat_gff.gff", sep = "\t")

    d_accepted = {}

    for interval_number in tqdm(df_ref.index.tolist(), desc = "Checking each interval in each genome"):
        d_accepted[interval_number] = {}
        for genome in list_genomes:
            if genome == ref_genome:
                d_accepted[interval_number][genome] = "Accepted"
            else :
                core1, core2 = get_core(df_ref.loc[interval_number, "left_core_family"], df_ref.loc[interval_number, "right_core_family"], genome, df_core)
                if check_contig(core1, core2) == False:
                    d_accepted[interval_number][genome] = "rejected (differents contigs)"
                else :
                    if check_rearrangement(core1, core2, df_gff, df_core, genome) == True:
                        d_accepted[interval_number][genome] = "Accepted"
                    else :
                        d_accepted[interval_number][genome] = "rejected (rearrangments)"


    df_accepted = pd.DataFrame.from_dict(d_accepted, orient = "index")
    list_columns = df_accepted.columns.tolist()
    df_accepted["n_accepted"] = df_accepted.apply(lambda x: count_accepted(x, list_columns), axis = 1)
    df_accepted = df_accepted.reset_index().rename(columns = {"index": "interval_number"})
    df_accepted.to_csv(f"{outdir}/1-core_intervals/all_genomes_check_corespot.tsv", sep = "\t", index = False)




def recreate_intervals(ref_genome, list_genomes, core_file, outdir):

    df_core = pd.read_csv(core_file, sep = "\t", names = ["gene_name", "core_family", "genome_name", "contig", "gene_number", "left_coordinate", "right_coordinate", "strand"], skiprows = 1)   
    df_gff = pd.read_csv(f"{outdir}/tmp/concat_gff.gff", sep = "\t")
    df_gff = df_gff[df_gff["type"] == "CDS"].reset_index(drop = True)
    df_check = pd.read_csv(f"{outdir}/1-core_intervals/all_genomes_check_corespot.tsv", sep = "\t")
    df_ref_intervals = pd.read_csv(f"{outdir}/1-core_intervals/Ref_{ref_genome}_intervals.tsv", sep = "\t")

    d_intervals = {}
    for interval_number in tqdm(df_ref_intervals["interval_number"].tolist(), desc = "Recreating all intervals in all genomes (if validated)"):
        d_intervals[interval_number] = {}
        for genome in list_genomes:
            if df_check.loc[interval_number,genome] == "Accepted" :
                core1, core2 = get_core(df_ref_intervals.loc[interval_number, "left_core_family"], df_ref_intervals.loc[interval_number, "right_core_family"], genome, df_core)
                d_intervals[interval_number][genome] = get_interval_genes(core1, core2, df_gff, genome)
            else :
                d_intervals[interval_number][genome] = "Rejected"
    
    df_intervals = pd.DataFrame.from_dict(d_intervals, orient = "index")
    df_intervals = df_intervals.reset_index().rename(columns = {"index": "interval_number"})
    df_intervals.to_csv(f"{outdir}/1-core_intervals/all_genomes_intervals.tsv", sep = "\t", index = False)


def get_core(core_family1, core_family2, genome, df_core):

    df_core2search = df_core[df_core["genome_name"] == genome]
    new_core1 = df_core2search[df_core2search["core_family"] == core_family1]["gene_name"].values[0]
    new_core2 = df_core2search[df_core2search["core_family"] == core_family2]["gene_name"].values[0]

    return new_core1, new_core2


def check_contig(core1, core2) :
    #return True if both core genes are on the same contig else False
    if core1.rsplit("_", 1)[0][:-1] == core2.rsplit("_", 1)[0][:-1]:
        return True
    else :
        return False


def check_rearrangement(core1, core2, df_gff, df_core, genome):

    # get the list of genes within the interval (core excluded)
    df_tmp = df_gff[df_gff["genome"] == genome].reset_index(drop = True)
    list_core = df_core[df_core["genome_name"] == genome]["gene_name"].tolist()
    list_core.remove(core1)
    list_core.remove(core2)

    index_core1 = df_tmp[df_tmp["id"] == core1].index[0]
    index_core2 = df_tmp[df_tmp["id"] == core2].index[0]
    
    list_genes_interval = df_tmp.iloc[min(index_core1, index_core2)+1:max(index_core1, index_core2)]["id"].tolist()

    #now checking if there are others core genes within the list_genes_interval
    for gene in list_genes_interval:
        if gene in list_core:
            return False
    
    return True


def count_accepted(row, list_columns):
    count = 0
    for columns in list_columns:
        if row[columns] == "Accepted":
            count += 1

    return count

def get_interval_genes(core1, core2, df_gff, genome):

    df_tmp = df_gff[df_gff["genome"] == genome].reset_index(drop = True)
    index_core1 = df_tmp[df_tmp["id"] == core1].index[0]
    index_core2 = df_tmp[df_tmp["id"] == core2].index[0]
    return (",").join(df_tmp.iloc[min(index_core1, index_core2)+1:max(index_core1, index_core2)]["id"].tolist())

