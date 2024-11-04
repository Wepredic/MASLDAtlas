def quality_control(adata,organism,sc,pd,np,mad, sns):
    
    # Calculating the qc metrics for the anndata
    sc.pp.calculate_qc_metrics(adata,
                                inplace=True, percent_top=[20], log1p=True)

    # Printing qc metrics before filtering
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], 
             jitter=0.4, multi_panel=True)

    # Filtering cells that has less than 200 counts
    sc.pp.filter_cells(adata, min_counts = 200)
    
    # Filtering genes that are not expressed at least 3 cells
    sc.pp.filter_genes(adata, min_cells = 3)

    # Adding annotation for mitochondrial genes 
    if organism == "zebrafish" or organism == "mouse":
        adata.var['mt'] = adata.var.index.str.startswith('mt-')
    elif organism == "human":
        adata.var['mt'] = adata.var.index.str.startswith('MT-')
    
    # Reading lists of ribosomal genes
    if organism == "zebrafish":
        ribo_genes = pd.read_csv("zebrafish_ribosomal.csv", header = None)
    elif organism == "human":
        ribo_genes = pd.read_csv("human_ribosomal.csv", header = None)
    elif organism == "mouse":
        ribo_genes = pd.read_csv("mouse_ribosomal.csv", header = None)

    # Adding annotation for ribosomal genes
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)

    # Calculating the qc metrics with mt- and ribosomal genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=True, inplace=True)
    
    # Median absolute deviations calculations
    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * mad(M)) | (
            np.median(M) + nmads * mad(M) < M
        )
        return outlier

    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    p2 = sc.pl.violin(adata, "pct_counts_mt")
    p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    # Filtering dead cells by mt pct > 20
    adata = adata[adata.obs.pct_counts_mt < 20]

    # 5 MAD outliers filtered out
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5)
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5))
    adata.obs.outlier.value_counts()

    # filtering 3 mads and mt pct exceeding 8 are filtered out
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 8)
    adata.obs.mt_outlier.value_counts()

    # Reporting # outlier cells
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    # Printing qc metrics after filtering
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)

    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    p2 = sc.pl.violin(adata, "pct_counts_mt")
    p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")


    return adata