def scrublet_scanpy(adata_path, out_path, batch_key = 'None', expected_doublet_rate=0.05):
  import scanpy as sc
  import pandas as pd
  adata = sc.read_h5ad(adata_path)
  if batch_key == 'None':
    batch_key = None
  adata = sc.external.pp.scrublet(adata, expected_doublet_rate = expected_doublet_rate, batch_key = batch_key, copy = True)
  df = adata.obs[['Cellname','doublet_score','predicted_doublet','Sample']]
  df.to_csv(out_path)
