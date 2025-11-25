import pandas as pd
import numpy as np
import scprep
import meld
import sklearn.mixture
import os
import argparse
def classes_Meld_py(metapath):
  # making sure plots & clusters are reproducible
  np.random.seed(42)
  metadata = pd.read_csv(metapath)
  print('Calculating GaussianMixture...')
  mixture_model = sklearn.mixture.GaussianMixture(n_components=3)
  likeh_array = metadata['stimulated_likelihood'].values.reshape(-1,1)
  print('fitting predict...')
  classes = mixture_model.fit_predict(likeh_array)
  classes = scprep.utils.sort_clusters_by_values(classes, metadata['stimulated_likelihood'])
  return(classes)
parser = argparse.ArgumentParser(
                    prog='Meld sub-function',
                    description='likehood to classes',
                    epilog='')
parser.add_argument('--inputpath')
parser.add_argument('--outputpath')
args = parser.parse_args()

classes = classes_Meld_py(args.inputpath)
DF = pd.DataFrame(classes)
DF.to_csv(args.outputpath)
