def meta_Meld_py(i, j, val, dim, cell_idx, feature_idx, group, test_method):
  import pandas as pd
  import numpy as np
  import graphtools as gt
  #import phate
  #import magic
  import scprep
  import meld
  import scipy
  #import seaborn as sns
  import scipy.io
  import os
  from scipy.sparse import csc_matrix
  # making sure plots & clusters are reproducible
  np.random.seed(42)
  data = csc_matrix((val, (i, j)), shape = dim)
  data.index = cell_idx
  data.columns = feature_idx
  #data = scprep.normalize.library_size_normalize(data)
  #data = scprep.transform.sqrt(data)
  metalist = [cell_idx, group]
  metadata = pd.DataFrame(metalist, index=['cell', 'condition']).T
  metadata.index = cell_idx
  meld_op = meld.MELD()
  sample_densities = meld_op.fit_transform(data, metadata['condition'])
  sample_densities.index = cell_idx
  print(data.shape)
  metadata['stimulated_likelihood'] = meld.utils.normalize_densities(sample_densities)['stimulated']
  return(metadata)
