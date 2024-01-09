# Load R libraries and source code
```
library(metacell)
source("xboc_functions.R")
```

# Run metacell clustering
Note: You can also start from the final metacell clustering object (provided in the "input_data/" folder, see next section).

```
##Initialize metacell database
scdb_init("mc_db",force_reinit=T)

##Add gene statistics
mcell_add_gene_stat(gstat_id="gstat", mat_id="mat_nobact_100", force=T)

##Load blacklisted genes
bl=scan("input_data/xboc_bl_v3",what="")

##Select variable markers for clustering
mcell_gset_filter_multi("gstat","clust_markers",T_tot=100,T_top3=2,T_szcor=-0.06,T_vm=2,T_niche=0.05,force_new=T,blacklist=bl)

##Plot selected marker statistics
mcell_plot_gstats(gstat_id="gstat", gset_id="clust_markers")

##Generate knn graph from matrix
mcell_add_cgraph_from_mat_bknn(mat_id="mat_nobact_100",gset_id="clust_markers",graph_id="graphk150",K=150,dsamp=F)

##Generate sc coclustering matrix by multiple metacell resamplings
mcell_add_cgraph_from_mat_bknn(mat_id="mat_nobact_100",gset_id="clust_markers",graph_id="graphk150",K=150,dsamp=F)

##Generate final metacell solution based on the coclustering matrix
mcell_mc_from_coclust_balanced("coc1000_k150_min30",mat="mat_nobact_100",mc_id="mc_K150_alpha2",K=30,min_mc_size=30,alpha=2)

##OPTIONAL: generate metacell 2D projection
library(tgconfig)
override_params("input_data/xboc_metacell_params.yaml","metacell")
mcell_mc2d_force_knn(mc2d_id="2dproj",mc_id="mc_K150_alpha2", graph_id="graphk150")
mcell_mc2d_plot(mc2d_id="2dproj")
mcell_mc2d_plot_by_factor("2dproj","mat_nobact_100","dataset",single_plot=T)

##Note: Bad metacells (low UMIs, no markers, etc) were filtered as described in the Methods section of the original article.
```

# Data visualization
```
##Initialize metacell database and load metacell and UMI matrix objects
scdb_init("mc_db",force_reinit=T)
mat=scdb_mat("mat_nobact_100")
mc=scdb_mc("mc_filt3")

##Load cell type annotation table and colorcode
xboc_ct_info=scr_load_cell_type_table("cell_type_annotation.txt",mc)
