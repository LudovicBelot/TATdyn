import pandas as pd
import sys
import os


#python script/sorted_TA_rows.py input/list_TA/TA_sorted_list.tsv.csv results/2023-01-26_TATdyn_6_full_photo_nopb/all_results_with_user_id_spots_only.tsv results/2023-01-26_TATdyn_6_full_photo_nopb/23-01-30_sorted_TA

def main():

    #small messy script to sort the results only for row with TAs (and sort them by experimental results)
    #Note here, we are considering the dynamism of a TA and the spot associated, if multiple TAs are located on the same spot, such spot will be represented multiple times in the results

    #first get the TA classified with the experimental results
    df_TA = pd.read_csv(sys.argv[1], sep = ";", comment = "#")
    df_dyn = pd.read_csv(sys.argv[2], sep = "\t", comment = "#")
    outdir = sys.argv[3]

    try :
        os.mkdir(outdir)
    except FileExistsError:
        print()

    d = {}
    for column in df_TA.columns.tolist():
        df_tmp = pd.DataFrame()
        list_elements = [x for x in df_TA[df_TA[column].notnull()][column].tolist() if x]
        for element in list_elements:
            if df_dyn[df_dyn["user_id"].str.contains(element )].empty == False:
                df_tmp = pd.concat([df_tmp, df_dyn[df_dyn["user_id"].str.contains(element )]])
    
        d[column] = df_tmp

                
    for k,v in d.items():
        v["Spot_lentgh"] = v.apply(lambda x: max(int(x["ref_left_c"]),int(x["ref_right_c"])) - min(int(x["ref_left_c"]),int(x["ref_right_c"])) , axis = 1) #Added the length of the spot
        
        v.to_csv(f"{outdir}/TAT_dyn_list_{k}.tsv", sep = "\t", index = False)


















if __name__ == "__main__":
    main()