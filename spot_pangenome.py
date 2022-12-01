import pandas as pd
import os
from tqdm import tqdm
from Bio import SeqIO

def mmseqs_align(outdir, list_genomes, n_threads):

    df = pd.read_csv(f"{outdir}/1-core_intervals/all_genomes_intervals.tsv", sep ="\t")
    #keeping only spots for which there are at least one gene
    df_intervals = df.loc[((df[list_genomes] != "Rejected") & (df[list_genomes].isna() == False)).any(axis = 1)]
    prt_seqio = SeqIO.index(f"{outdir}/tmp/cat_proteins.prt", "fasta")

    d_interval_genes = {} # k = interval_number, v = {gene_family_number: [list_of_genes_in_all_genomes_clusterized]}
    #Now for every spot, we get all proteins seqs in one file and run MMseqs2 (if n_genes >1)
    for interval in tqdm(df_intervals.iterrows(), total = df_intervals.shape[0], desc = "Running MMseqs clustering for every intervals"):
        list_genes2cluster = row2list(interval,list_genomes)
        if len(list_genes2cluster) > 1:
            d_interval_genes[interval[1]["interval_number"]] = interval_cluster(list_genes2cluster, prt_seqio, outdir, n_threads)
        elif len(list_genes2cluster) == 1:
            d_interval_genes[interval[1]["interval_number"]] = {1: list_genes2cluster[0]}
    
    d_interval_genes, d2_save = format_dic(d_interval_genes)

    #saving the file for manual check if needed
    df2_save = pd.DataFrame.from_dict(d2_save, orient = "index")
    df2_save.to_csv(f"{outdir}/2-spot_pangenome/spot_pangenome_clusterized.tsv", sep = "\t", index = False)

    #Now creating a tsv file with row = {Interval_number}_{gene_family_number}, columns = one column per genome
    # The value represents the number of gene homologs part of each family in each genome 
    # Note in case we rejected the genome before: the value will be "?" for the studied genome
    df_table = create_cluster_table(d_interval_genes, list_genomes, df)
    df_table.to_csv(f"{outdir}//2-spot_pangenome/spot_pangenome_summary.table", sep = "\t")





def row2list(row, list_genomes):
    #return a list of all genes within the same interval in all genomes
    res = []
    for genome in list_genomes:
        if row[1][genome] != None and row[1][genome] != "Rejected":
            res += row[1][genome].split(",")
    
    return res


def interval_cluster(list_genes2cluster, prt_seqio, outdir, n_threads):

    #return a dictionnary with for the given interval : k = family_gene_number, v = [list_of_genes_in_all_genomes_clusterized]
    #first writing all seqs from the same interval into a tmp file to run MMseqs2
    str_interval = ""
    for gene in list_genes2cluster:
        str_interval += f">{gene}\n{prt_seqio[gene].seq}\n"
    
    with open(f"{outdir}/tmp/tmp_spot2cluster.fasta", "w") as f:
        f.write(str_interval)
    
    #running MMseqs2
    os.system(f"mmseqs easy-cluster {outdir}/tmp/tmp_spot2cluster.fasta {outdir}/tmp/spot_cluster {outdir}/tmp/mmseqs --min-seq-id 0.8 -v 0 --threads {n_threads} --remove-tmp-files")
    #Sorting every gene family per spot
    df_cluster = pd.read_csv(f"{outdir}/tmp/spot_cluster_cluster.tsv", sep = "\t", names = ["ref", "target"])
    n_family = 1
    d_families = {}
    for ref_gene in df_cluster["ref"].drop_duplicates():
        d_families[n_family] = df_cluster[df_cluster["ref"] == ref_gene]["target"].tolist()
        n_family +=1
    
    return d_families


def format_dic(d2format):

    #small function which return two dict as follow
    # one useful for next operations: k = {interval_number}_{gene_family_number}, v = [list_of_genes_in_all_genomes_clusterized]
    #and one just to save an intermediate file
    d = {}
    d2save = {}
    n = 0
    for interval_number, d_spot in d2format.items():
        for family_number, list_genes in d_spot.items():
            d[f"{interval_number}_{family_number}"] = list_genes
            d2save[n] = {"Spotnumber_family":f"{interval_number}_{family_number}", "list_genes":(",").join(list_genes)}
            n+=1
    
    return d, d2save

def create_cluster_table(d_interval_genes, list_genomes, df):
    #return a df with row = {Interval_number}_{gene_family_number}, columns = one column per genome
    # The value represents the number of gene homologs part of each family in each genome 
    # Note in case we rejected the genome before: the value will be "?" for the studied genome

    d_res = {}
    n_index = 0
    for gene_family, list_genes in d_interval_genes.items():
        n_index += 1
        d_res[gene_family] = {}
        for genome in list_genomes:
            if df.loc[int(gene_family.split("_")[0]), genome] != "Rejected":
                d_res[gene_family][genome] = count_gene_in_family(list_genes, genome)
            elif df.loc[int(gene_family.split("_")[0]), genome] == "Rejected":
                d_res[gene_family][genome] = "?"

    return pd.DataFrame.from_dict(d_res, orient = "index")



def count_gene_in_family(list_genes, genome):
    #small function which count the number of genes part of the given genome in the list_genes
    
    count = 0 
    for gene in list_genes:
        if gene.rsplit(".",1)[0] == genome:
            count += 1
    
    return count 