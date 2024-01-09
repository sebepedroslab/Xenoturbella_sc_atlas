# libraries
import pandas as pd
import numpy as np
import sys
import samalg
import scanpy as sc

# output
out_fn = "results_alignment_samap/"

## Preprocess ##

# # reference species
# spi = "Xboc"

# # read
# print("# loading %s" % (spi))
# mat_i = sc.read_mtx("data/expression/%s_UMI_table.mtx" % spi)
# sct_i = pd.read_csv("data/expression/%s_cellID_cell_type" % spi, sep="\t", header = None, names = ["cell","celltype"])
# gen_i = pd.read_csv("data/expression/%s_features.tsv" % spi, sep="\t", header = None, names = ["gene"])

# # create SAM objects for concatenated dataset
# print("# run SAM %s" % (spi))
# sam_i = samalg.SAM(counts=[mat_i.X.transpose(), gen_i["gene"].values,sct_i["cell"].values ])
# sam_i.preprocess_data()
# for i in range(sct_i.shape[1]):
# 	sam_i.adata.obs[sct_i.columns[i]] = sct_i.iloc[:,i].values
# sam_i.run()

# # save precomputed sam object
# sam_i.save_anndata("%s/data.%s.sam_object.h5ad" % (out_fn, spi))


# query species
for spj in ["Xboc","Hmia","Spur","Isopu","Nvec"]:
	
	# read
	print("# loading %s" % (spj))
	mat_j = sc.read_mtx("data/expression/%s_UMI_table.mtx" % spj)
	sct_j = pd.read_csv("data/expression/%s_cellID_cell_type" % spj, sep="\t", header = None, names = ["cell","celltype"])
	gen_j = pd.read_csv("data/expression/%s_features.tsv" % spj, sep="\t", names = ["gene"])

	# create SAM objects for concatenated dataset
	print("# run SAM %s" % (spj))
	sam_j = samalg.SAM(counts=[mat_j.X.transpose(), gen_j["gene"].values,sct_j["cell"].values ])
	sam_j.preprocess_data()
	for i in range(sct_j.shape[1]):
		sam_j.adata.obs[sct_j.columns[i]] = sct_j.iloc[:,i].values
	sam_j.run()

	# save precomputed sam object
	sam_j.save_anndata("%s/data.%s.sam_object.h5ad" % (out_fn, spj))

print("all done!")