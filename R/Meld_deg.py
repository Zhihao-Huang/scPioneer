def deg_Meld_py(i, j, val, dim, cell_idx, feature_idx, classes, test_method):
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
  data = pd.DataFrame(data.toarray())
  data.index = cell_idx
  data.columns = feature_idx
  classes = np.asarray(classes)
  # DEG
  # Here we take the cells from the depleted and enriched groups (class 0 and 2, respectively)
  # and use them as input to the differential expresssion test
  de_results = de.test.two_sample(data[np.isin(classes, [0,2])].values, 
                                grouping=classes[np.isin(classes, [0,2])],
                                gene_names=data.columns,
                                test=test_method,
                                noise_model='nb').summary()
  return(de_results)
