# libraries
import sys
import numpy as np
import pandas as pd
import samap
import scanpy as sc
import samalg
from samap.mapping import SAMAP
from samap.utils import save_samap, load_samap

# output
out_fn = "results_alignment_samap/"

# list of pairwise species comparisons
com_list = [ ["Xboc","Spur"], ["Xboc","Hmia"], ["Xboc","Isopu"], ["Xboc","Nvec"] ]
com_list = [ ["Xboc","Nvec"] ]

for n,com in enumerate(com_list):
	
	# query species data
	sid1 = com[0]
	sid2 = com[1]
	
	### Load data ###
	
	# reload SAM objects for species 1
	print("# loading %s data" % sid1)
	sam1 = samalg.SAM()
	sam1.load_data("%s/data.%s.sam_object.h5ad" % (out_fn, sid1))
	sc.pp.highly_variable_genes(sam1.adata)
	# sam1.leiden_clustering(res=3)

	# reload SAM objects for species 2
	print("# loading %s data" % sid2)
	sam2 = samalg.SAM()
	sam2.load_data("%s/data.%s.sam_object.h5ad" % (out_fn, sid2))
	sc.pp.highly_variable_genes(sam2.adata)
	# sam2.leiden_clustering(res=3)

	# run samap, all homlogs with leiden clusters
	print("# run SAMAP %s-%s" % (sid1,sid2))
	samm = SAMAP({ sid1:sam1, sid2:sam2 }, f_maps = "data/blast/", keys = {sid1:"celltype", sid2:"celltype"})
	samm.run()
	save_samap(samm, "%s/data.pair.%s-%s" % (out_fn, sid1, sid2))

print("all done!")
