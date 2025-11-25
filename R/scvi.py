def DEG_scvi_py(i, j, val, dim, cell_idx, feature_idx, meta, batch_key, n_cores):
  import os
  import tempfile
  import matplotlib.pyplot as plt
  import numpy as np
  import pandas as pd
  #import plotnine as p9
  #import scanpy as sc
  import scvi
  #import seaborn as sns
  import torch
  import anndata as ad
  from scipy.sparse import csc_matrix
  np.random.seed(42)
  data = csc_matrix((val, (i, j)), shape = dim)
  if batch_key == 'None':
    batch_key = None
  adata = ad.AnnData(data)
  adata.obs_names = cell_idx
  adata.var_names = feature_idx
  print(adata)
  for i in meta.columns:
    adata.obs[i] = meta[i]
  scvi.data.poisson_gene_selection(adata)
  print(adata.var.head())
  adata = adata[:, adata.var["highly_variable"]]
  adata.layers["counts"] = adata.X.copy().tocsr()  # converts to CSR format, preserve counts
  scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key = batch_key)  # prepare data for scVI
  ## Define and train the model
  model = scvi.model.SCVI(
    adata, gene_likelihood="nb"
  )  # We use Negative Binomial count likelihoods, following Boyeau et al., 2023.
  model.train(
    accelerator="cpu", devices='auto', #devices=int(n_cores),
    check_val_every_n_epoch=1,
    max_epochs=400,
    early_stopping=True,
    early_stopping_patience=20,
    early_stopping_monitor="elbo_validation")
  condition_1 = "stimulated"
  cell_idx1 = adata.obs["condition"] == condition_1
  print(sum(cell_idx1), "condition", condition_1)
  condition_2 = "unstimulated"
  cell_idx2 = adata.obs["condition"] == condition_2
  print(sum(cell_idx2), "condition", condition_2)
  de_change = model.differential_expression(idx1=cell_idx1, idx2=cell_idx2,weights="uniform", batch_correction=True)
  return(de_change)
