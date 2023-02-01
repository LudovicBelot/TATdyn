import argparse
from datetime import date, datetime
import os.path

import utilities
import recreate_interval
import spot_pangenome
import HTgenes
import HCspot
import Spot_dissimilarity
import TAT_defense_finder

#commandline example : 
# python script/main.py -d input/full_photo -c input/full_photo/fully_assembled_photo_corefeatures.tsv -o results -t input/full_photo/tree_full_photo_rooted_asymbiotica.nwk -r Ph22.1222.00005 --threads 0 --mygenes input/TT01_TA_coordinates.tsv
# python script/main.py -c input/core_genome_photorhabdus_genus_features.csv --indir input/panacota_photo --outdir results --tree input/panacota_photo/photo_genus.treefile --mygenes input/list_TA/TT01_TA_coordinates.tsv --threads 0


def main():

    args = get_args()
    today_date = date.today()
    utilities.check_dir(args.outdir, today_date)
    ref_genome, list_genomes = utilities.ref_genome(args.indir, args.ref)
    utilities.combine_input(args.indir, f"{args.outdir}/{today_date}_TATdyn", list_genomes)
    n_threads = utilities.threads2use(args.threads)
    n_cores = utilities.core2use(args.cpu)
    #adding names for each non-terminals nodes if they do not exist, useful later to determine which genes are HTgenes
    utilities.name_tree(args.tree, f"{args.outdir}/{today_date}_TATdyn")

    #First we need to recreate each interval in each genome
    #Note: The first genome in your folder will be considered as reference or you can use -r --ref {genome_name} to set the given genome as the reference
    #This is important because the order of the interval will be based on this genome
    recreate_interval.ref_interval(ref_genome, args.core, f"{args.outdir}/{today_date}_TATdyn")
    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/1-core_intervals/all_genomes_check_corespot.tsv") == False:
        recreate_interval.check_interval(ref_genome, list_genomes, args.core, f"{args.outdir}/{today_date}_TATdyn")
    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/1-core_intervals/all_genomes_intervals.tsv") == False:
        recreate_interval.recreate_intervals(ref_genome, list_genomes, args.core, f"{args.outdir}/{today_date}_TATdyn")


    #Using MMseqs2 to cluster all proteomes together useful to create the Count model rates (birth death model)
    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/2-spot_pangenome/whole_pangenome_summary.table") == False:
        spot_pangenome.mmseqs_whole_align(f"{args.outdir}/{today_date}_TATdyn", list_genomes, n_threads)


    #Now using MMseqs2, we align all proteins in each pangenome spot
    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/2-spot_pangenome/spot_pangenome_summary.table") == False:
        spot_pangenome.mmseqs_align(f"{args.outdir}/{today_date}_TATdyn", list_genomes, n_threads)

    #running Count to determine which genes are Horizontally transfered
    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/2-spot_pangenome/Final_spot_pangenome_HTg.tsv") == False:
        HTgenes.run_Count(f"{args.outdir}/{today_date}_TATdyn")
        HTgenes.determine_HTgenes(f"{args.outdir}/{today_date}_TATdyn")
    
    n_HTevents, m_spots = HTgenes.HTg_per_spot(f"{args.outdir}/{today_date}_TATdyn")


    #running a simulations of random distribution of the n HTevents within the m Spots
    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/3-HC_analysis/HC_summary.txt") == False:
        t95 = HCspot.HC_sim(n_HTevents, m_spots, f"{args.outdir}/{today_date}_TATdyn")
        HCspot.HC_spot(t95, f"{args.outdir}/{today_date}_TATdyn")

    # Now determining the dissimilarity index 
    # Ref = "Partioning the turnover and nestedness components of beta diversity, Andrés Baselga, 2010 Global Ecology and Biogeography"
    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/4-Spot_dissimilarity/Final_spot_dissimilarity_index.tsv") == False:
        list_combinaisons, index_genome = Spot_dissimilarity.prepare_data(f"{args.outdir}/{today_date}_TATdyn", n_core = n_cores)
        Spot_dissimilarity.calculate_dissimilarity(list_combinaisons, index_genome, f"{args.outdir}/{today_date}_TATdyn")

    if os.path.exists(f"{args.outdir}/{today_date}_TATdyn/5-optional_defense-finder/1_spot_pangenome_defense_system.tsv") == False and args.defense_finder == True:
        os.system(f"cat {args.indir}/gff3/*.gff > {args.outdir}/{today_date}_TATdyn/tmp/cat4defense_gff")
        TAT_defense_finder.spot_defense_finder(f"{args.outdir}/{today_date}_TATdyn/1-core_intervals/all_genomes_intervals.tsv", 
                                                f"{args.outdir}/{today_date}_TATdyn/tmp/cat_proteins.prt", 
                                                f"{args.outdir}/{today_date}_TATdyn/tmp/cat4defense_gff", 
                                                f"{args.outdir}/{today_date}_TATdyn/5-optional_defense-finder", 
                                                n_core=n_cores)

    #combine the results in one final file
    utilities.combine_results(f"{args.outdir}/{today_date}_TATdyn", user_gene = args.mygenes, TAT_defense = args.defense_finder)


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("--indir", "-d",
                        help = " (REQUIRED) Directory in which are located the proteins/gff/replicons generated by PanACoTA annotate", required = True)
    parser.add_argument("--core", "-c",
                help = "(REQUIRED) Core features files generated using PanACoTA corepers (strict 100%%) followed by TAT_core2df.py script",
                required = True)
    parser.add_argument("--outdir", "-o",
                help ="(REQUIRED) Out directory in which the results will be stored",
                required = True)
    parser.add_argument("--tree", "-t",
                        help = "(REQUIRED) Treefile generated by PanACoTA tree module (or software such as Iqtree,...) in newick format. Needed if the non-terminals nodes are not named",
                        required = True
                        )
    parser.add_argument("-r", "--ref",
                        help = "If you want to specify which genome in you dataset will serve as reference (default: the first in your dataset), \n Note: At least this genome needs to be fully assembled",
                        default = None
                        )
    parser.add_argument("--mygenes", "-g",
                        help = "Tsv file with the genes you want",
                        default= None
                        )
    parser.add_argument("--defense_finder",
                        help = "If you want to also research potential antiviral system in every spot, it will use Defense-finder developed by F. Tesson et al, 2022",
                        action='store_true',
                        default= False
                        )
    parser.add_argument("--cpu",
                        help ="Number of cores to use (default = 1), 0 for the max core to use",
                        default = 1
                        )
    parser.add_argument("--threads",
                        help = "Number of threads to use (default = 1) 0 for the max threads to use",
                        default = 1
                        )

    args = parser.parse_args()
    return args






if __name__ == "__main__":
    main()