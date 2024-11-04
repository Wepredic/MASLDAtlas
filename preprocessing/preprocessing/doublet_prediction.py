def doublet_prediction(adata,sc,scvi):

    # Copying the anndata since it will be subsetted for highly
    # variable genes
    doublet_removal = adata.copy()


    # Calculation of Highly Variable Genes
    sc.pp.highly_variable_genes(doublet_removal, n_top_genes = 2000,
                                 subset = True, flavor = 'seurat_v3')

    # Setting up & training SCVI model
    scvi.model.SCVI.setup_anndata(doublet_removal)
    vae = scvi.model.SCVI(doublet_removal)
    vae.train()

    # Setting up & training SOLO model    
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()

    # Getting prediction scores for singlets and doublets,
    # Getting indexes of cells which have >1 difference between singlet
    # and doublet prediction scores
    df = solo.predict()
    df['prediction'] = solo.predict(soft = False)
    df.index = df.index.map(lambda x:x)
    df.groupby('prediction').count()
    df['difference'] = df.doublet - df.singlet
    doublets = df[(df.prediction == 'doublet') & (df.difference > 1)]




    return doublets