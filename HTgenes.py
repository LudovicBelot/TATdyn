import pandas as pd





def determine_HTgenes(countfile, outdir):
    #function which determine which gene from a given pangenome is considered as Horizontal Transfer genes (HTg)
    # This using count and the birth - loss - duplication model

    #First get a list of the columns names in the file
    n = 0
    with open(countfile, "r") as f:
        for line in f:
            n+=1
            if n == 3:
                columns_names = [x.strip() for x in line.split("\t") if x]
                break


    df_count = pd.read_csv(countfile, sep = "\t", comment = "#", names = columns_names)
    print(df_count[df_count["PhAl.1022.00028"] != 0]["# Family"].drop_duplicates().tolist())
    print(len(df_count[df_count["PhAl.1022.00028"] != 0]["# Family"].drop_duplicates().tolist()))
    """
    list_columns2keep = [x for x in columns_names if ":gain" in x]
    list_columns2keep.insert(0,"# Family")
    df_count = df_count[list_columns2keep]
    list_columns2keep.remove("# Family")
    print(df_count[list_columns2keep])
    #df_count["n_acquisitions"] = df_count[df_count[list_columns2keep] > 0.95].count()
    print(df_count[df_count[list_columns2keep] > 0.95].count())
    print(df_count["PhAL.1022.00028"])
    """