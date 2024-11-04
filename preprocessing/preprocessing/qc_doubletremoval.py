def qc_doubletremoval(folder_path,file_type, organism,sc,scvi,pd,np, mad, sns):
    from quality_control import quality_control
    from doublet_prediction import doublet_prediction

    print(folder_path)


    # The rawfiles mostly mtx, but in any case I added them,
    # there are more file types, check scanpy for more read functions
    if file_type == "mtx":
        adata = sc.read_10x_mtx(folder_path)
    elif file_type == "h5ad":
        adata = sc.read_h5ad(folder_path)
    elif file_type == "csv":
        adata = sc.read_csv(folder_path)
    elif file_type == "excel":
        adata = sc.read_excel(folder_path)
    elif file_type == "h5":
        adata = sc.read_h5(folder_path)

    # filtering qc
    adata = quality_control(adata,organism,sc,pd,np, mad, sns)
    
    # getting indexes for doublet predicted cells
    doublets = doublet_prediction(adata,sc,scvi)

    # removing doublets 
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]
    
    #final qc metrics of filtered data
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)    
    
    
    return adata