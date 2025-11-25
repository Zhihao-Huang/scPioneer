def all_Meld_py(i, j, val, dim, cell_idx, feature_idx, group, test_method):
  import pandas as pd
  import numpy as np
  import graphtools as gt
  #import phate
  #import magic
  import scprep
  import meld
  import scipy
  import seaborn as sns
  import diffxpy.api as de
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
  print(metadata.head())
  metadata['stimulated_likelihood'] = meld.utils.normalize_densities(sample_densities)['stimulated']
  print('Calculating GaussianMixture...')
  import sklearn.mixture
  mixture_model = sklearn.mixture.GaussianMixture(n_components=3)
  likeh_array = metadata['stimulated_likelihood'].values.reshape(-1,1)
  print('fitting predict...')
  classes = mixture_model.fit_predict(likeh_array)
  classes = scprep.utils.sort_clusters_by_values(classes, metadata['stimulated_likelihood'])
  if all(np.isin([0,2],np.unique(classes))):
    # DEG
    # Here we take the cells from the depleted and enriched groups (class 0 and 2, respectively)
    # and use them as input to the differential expresssion test
    data = pd.DataFrame(data.toarray())
    de_results = de.test.two_sample(data[np.isin(classes, [0,2])].values,
      grouping=classes[np.isin(classes, [0,2])],
      gene_names=data.columns,
      test=test_method).summary()
  else:
    print('Depleted or enriched group was absent in classes. Return empty.')
    de_results = ''
  return(classes, de_results)
