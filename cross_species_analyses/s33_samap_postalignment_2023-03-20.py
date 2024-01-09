# libraries
import sys
import numpy as np
import pandas as pd
import samap
from samap.utils import save_samap, load_samap, df_to_dict, substr
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from samalg import SAM
import scipy as sp

# output
out_fn = "results_alignment_samap/"

# list of pairwise species comparisons
com_list = [ ["Xboc","Isopu"], ["Xboc","Spur"], ["Xboc","Hmia"], ["Xboc","Nvec"] ]
com_list = [ ["Xboc","Nvec"] ]


# functions
def q(x):
	return np.array(list(x))

def smod_get_mapping_scores_top_f(sm, keys, f_top = 0):
	"""Calculate mapping scores
	Parameters
	----------
	sm: SAMAP object

	keys: dict, annotation vector keys for at least two species with species identifiers as the keys
		e.g. {'pl':'tissue','sc':'tissue'}

	n_top: int, optional, default 0
		If `n_top` is 0, average the alignment scores for all cells in a pair of clusters.
		Otherwise, average the alignment scores of the top `n_top` cells in a pair of clusters.
		Set this to non-zero if you suspect there to be subpopulations of your cell types mapping
		to distinct cell types in the other species.
	Returns
	-------
	D - table of highest mapping scores for cell types 
	A - pairwise table of mapping scores between cell types across species
	"""


	if len(list(keys.keys()))<len(list(sm.sams.keys())):
		samap = SAM(counts = sm.samap.adata[np.in1d(sm.samap.adata.obs['species'],list(keys.keys()))])
	else:
		samap=sm.samap

	clusters = []
	ix = np.unique(samap.adata.obs['species'],return_index=True)[1]
	skeys = q(samap.adata.obs['species'])[np.sort(ix)]

	for sid in skeys:
		clusters.append(q([sid+'_'+str(x) for x in sm.sams[sid].adata.obs[keys[sid]]]))

	cl = np.concatenate(clusters)
	l = "{}_mapping_scores".format(';'.join([keys[sid] for sid in skeys]))
	samap.adata.obs[l] = pd.Categorical(cl)

	CSIMth, clu = smod_compute_csim_top_f(samap, l, f_top = f_top, prepend = False)

	A = pd.DataFrame(data=CSIMth, index=clu, columns=clu)
	i = np.argsort(-A.values.max(0).flatten())
	H = []
	C = []
	for I in range(A.shape[1]):
		x = A.iloc[:, i[I]].sort_values(ascending=False)
		H.append(np.vstack((x.index, x.values)).T)
		C.append(A.columns[i[I]])
		C.append(A.columns[i[I]])
	H = np.hstack(H)
	D = pd.DataFrame(data=H, columns=[C, ["Cluster","Alignment score"]*(H.shape[1]//2)])
	return D, A

def smod_compute_csim_top_f(samap, key, X=None, prepend=True, f_top = 0):
	splabels = q(samap.adata.obs['species'])
	skeys = splabels[np.sort(np.unique(splabels,return_index=True)[1])]

	cl = []
	clu = []
	for sid in skeys:
		if prepend:
			cl.append(sid+'_'+q(samap.adata.obs[key])[samap.adata.obs['species']==sid].astype('str').astype('object'))
		else:
			cl.append(q(samap.adata.obs[key])[samap.adata.obs['species']==sid])            
		clu.append(np.unique(cl[-1]))

	clu = np.concatenate(clu)
	cl = np.concatenate(cl)

	CSIM = np.zeros((clu.size, clu.size))
	if X is None:
		X = samap.adata.obsp["connectivities"].copy()

	xi,yi = X.nonzero()
	spxi = splabels[xi]
	spyi = splabels[yi]

	filt = spxi!=spyi
	di = X.data[filt]
	xi = xi[filt]
	yi = yi[filt]

	px,py = xi,cl[yi]
	p = px.astype('str').astype('object')+';'+py.astype('object')

	A = pd.DataFrame(data=np.vstack((p, di)).T, columns=["x", "y"])
	valdict = df_to_dict(A, key_key="x", val_key="y")   
	cell_scores = [valdict[k].sum() for k in valdict.keys()]
	ixer = pd.Series(data=np.arange(clu.size),index=clu)
	if len(valdict.keys())>0:
		xc,yc = substr(list(valdict.keys()),';')
		xc = xc.astype('int')
		yc=ixer[yc].values
		cell_cluster_scores = sp.sparse.coo_matrix((cell_scores,(xc,yc)),shape=(X.shape[0],clu.size)).A

		for i, c in enumerate(clu):
			if f_top > 0:
				n_top = int(len(np.sort(cell_cluster_scores[cl==c],axis=0)) * f_top)
				CSIM[i, :] = np.sort(cell_cluster_scores[cl==c],axis=0)[-n_top:].mean(0)
			else:
				CSIM[i, :] = cell_cluster_scores[cl==c].mean(0)

		CSIM = np.stack((CSIM,CSIM.T),axis=2).max(2)
		CSIMth = CSIM / samap.adata.uns['mapping_K']
		return CSIMth,clu
	else:
		return np.zeros((clu.size, clu.size)), clu




# loop
for n,com in enumerate(com_list):
	
	# query species data
	sid1 = com[0]
	sid2 = com[1]
	
	### Load data ###
	
	# load samap
	print("%s-%s | SAMAP load object" % (sid1,sid2))
	samm = load_samap("%s/data.pair.%s-%s" % (out_fn, sid1, sid2))



	### Mapping scores ###
	
	# broad cell types
	print("%s-%s | SAMAP alignment scores, cell types" % (sid1,sid2))
	map_topct,map_score = smod_get_mapping_scores_top_f(samm, keys = {sid1:"celltype", sid2:"celltype"}, f_top = 0.5)
	map_score.to_csv("%s/samap.%s-%s.scores.cts.top100q.csv" % (out_fn, sid1, sid2), sep = "\t")
	map_topct.to_csv("%s/samap.%s-%s.topcts.cts.top100q.csv" % (out_fn, sid1, sid2), sep = "\t")
	
	print("%s-%s | SAMAP alignment scores, gene pairs" % (sid1,sid2))
	gpf = samap.analysis.GenePairFinder(samm, keys = {sid1:"celltype", sid2:"celltype"})
	gene_pairs = gpf.find_all(align_thr = 0.15)
	gene_pairs.to_csv("%s/samap.%s-%s.scores.cts.shared_genes.csv" % (out_fn, sid1, sid2), index = False)
		
	print("%s-%s | SAMAP alignment scores, plots" % (sid1,sid2))
	with PdfPages("%s/samap.%s-%s.scores.cts.plots.pdf" % (out_fn, sid1, sid2)) as pdf:
		samm.scatter(COLORS = { sid1:"blue", sid2:"gray" })
		pdf.savefig()
		dd = sc.pl.umap(
			samm.samap.adata,
			show=False,
			frameon = True,
			color='%s_celltype' % sid1,
			add_outline=False,
			legend_loc='on data',
			legend_fontsize=7,
			legend_fontoutline=2,
			title=sid1,
			palette='viridis')
		pdf.savefig()
		dd = sc.pl.umap(
			samm.samap.adata,
			show=False,
			frameon = True,
			color='%s_celltype' % sid2,
			add_outline=False,
			legend_loc='on data',
			legend_fontsize=7,
			legend_fontoutline=2,
			title=sid2,
			palette='viridis')
		pdf.savefig()
		dd = sc.pl.umap(
			samm.samap.adata,
			show=False,
			frameon = True,
			color='celltype;celltype_mapping_scores',
			add_outline=False,
			legend_loc='on data',
			legend_fontsize=7,
			legend_fontoutline=2,
			title=sid2,
			palette='viridis')
		pdf.savefig()
		plt.close()


print("all done!")
