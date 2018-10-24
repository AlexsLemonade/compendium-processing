import argparse
import pandas as pd
import numpy as np
import random
from fancyimpute import KNN, BiScaler, SoftImpute, IterativeSVD
from sklearn import preprocessing

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help = "Masked .pcl file")
args = parser.parse_args()
input_file = args.file 

# output file "base"
outfile = input_file.replace(".pcl", "").replace("data/masked", "imputed")

# set random seed 2 ways cause I'm not sure what's appropriate, my suspicion
# is numpy
random.seed(123)
np.random.seed(123)

# read in data and transpose
data = pd.read_csv(input_file, sep='\t', header=0, index_col=0, 
				   error_bad_lines=False)
new_data = data.copy()
transposed = new_data.T

# we'll need a matrix specifically for the biscaler transform, for SoftImpute
print("SoftImpute...")
transposed_mat = transposed.as_matrix()
biscaler = BiScaler()

# perform the scaling appropriate for this imputation strategy
transposed_normalized = biscaler.fit_transform(transposed_mat)

# the imputation itself
imputed_softimpute = SoftImpute().fit_transform(transposed_normalized)

# we don't want the transformed values and we want samples to be columns
inverse_softimpute = biscaler.inverse_transform(imputed_softimpute)
untransposed_softimpute = inverse_softimpute.transpose()

# prepare to write to file, back to DataFrame, return indices
softimpute_df = pd.DataFrame(untransposed_softimpute)
softimpute_df.index = data.index
softimpute_df.columns = data.columns.values

# write to a tab separated values file, but we'll use the .pcl file extension
softimpute_outfile = outfile + "_softimpute.pcl"
softimpute_df.to_csv(softimpute_outfile, sep='\t')

print("KNN...")
# KNN assumes that standard scaling has been performed
scaler = preprocessing.StandardScaler(copy=True)
scaler.fit(transposed)
scaled = pd.DataFrame(scaler.transform(transposed),
                      index=transposed.index,
                      columns=transposed.columns
)

# perform the imputation, setting k=10 as is standard for gene expression data
imputed_knn = KNN(k=10).fit_transform(scaled)

# inverse transformation -- we don't want the standard scores
inverse_knn = scaler.inverse_transform(imputed_knn)

# columns are samples
untransposed_knn = inverse_knn.transpose()

# write to file
knn_df = pd.DataFrame(untransposed_knn)
knn_df.index = data.index
knn_df.columns = data.columns.values
# not to be confused with the Sleipnir KNNImputer output
knn_outfile = outfile + "_KNN_fancyimpute.pcl"
knn_df.to_csv(knn_outfile, sep='\t')

print("IterativeSVD...")
# standard scaled
imputed_svd = IterativeSVD(rank=10).fit_transform(transposed)

# inverse transform
inverse_svd = scaler.inverse_transform(imputed_svd)

# columns are samples
untransposed_svd = imputed_svd.transpose()

# write to file
svd_df = pd.DataFrame(untransposed_svd)
svd_df.index = data.index
svd_df.columns = data.columns.values
# not to be confused with the Sleipnir KNNImputer output
svd_outfile = outfile + "_IterativeSVD.pcl"
svd_df.to_csv(svd_outfile, sep='\t')