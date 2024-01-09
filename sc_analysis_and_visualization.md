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

# Metacell data visualization
```
##Initialize metacell database and load metacell and UMI matrix objects
scdb_init("mc_db",force_reinit=T)
mat=scdb_mat("mat_nobact_100")
mc=scdb_mc("mc_filt3")

##Load cell type annotation table and colorcode
xboc_ct_info=scr_load_cell_type_table("input_data/cell_type_annotation.txt",mc)

##Generate metacell-level and cell type-level raw and normalized (umifrac) UMI matrices
mc_counts=scr_mc_gene_counts(mc,mat,5)
mc_umifrac=scr_mc_gene_umifrac(mc,mat,mc_counts,5)
ct_mc=scr_cell_type_fp("input_data/cell_type_annotation.txt",mc,mat)
ct_counts=scr_mc_gene_counts(ct_mc,mat,5)
ct_umifrac=scr_mc_gene_umifrac(ct_mc,mat,ct_counts,5)

##Figure 1C
mcell_mc2d_plot(mc2d_id="2dproj_v3")

##Figure 1D
scr_plot_cmod_markers_ct_colors(mc,mat,clust_ord=niche_order,output_file="Global_gene_expression.png",gene_annot_file="input_data/xboc_gene_annotation",per_clust_genes=30,height=10000,gene_min_fold=2,sn_table="input_data/cell_type_annotation.txt")

##Dotmap example (Figures 1H, 2B, 2F, 3A, 3E and S3A)
scr_dot_plot_map(mc_umifrac,mc,"input_data/Specific_lists/example_markers",xboc_ct_info,out_fn="example_dotmap")

##Figure S2 - TF heatmap
scr_barplot_heatmap_markers(mc,mat,mc_counts,"input_data/xboc_tfs",heatmap_file="TFs_heatmap.pdf",min_gene_fc=1.8,print_barplots=F,pmin=3.5,mc_col=xboc_ct_info$mc_color,sort_markers=T)
```

# Bacterial symbionts analysis
```
##Initialize metacell database and oad UMI matrix containing bacterial counts
scdb_init("mc_db",force_reinit=T)
mat_bact=scdb_mat("mat")
mat_no_bact=scdb_mat("mat_nobact_100")
mc=scdb_mc("mc_filt3")

#Count bacterial signal per metacell
chlamydia=colSums(as.matrix(mat_bact@mat[grepl("Chlamydia",rownames(mat_bact@mat)),names(mc@mc)]))
proteobact=colSums(as.matrix(mat_bact@mat[grepl("Proteobact",rownames(mat_bact@mat)),names(mc@mc)]))

bact_tot=colSums(rbind(chlamydia,proteobact))
chlamydia_tot=tapply(chlamydia[names(mc@mc)],mc@mc,sum)
proteobact_tot=tapply(proteobact[names(mc@mc)],mc@mc,sum)

mc_counts=scr_mc_gene_counts(mc,mat_no_bact,5)
mc_sizes=colSums(mc_counts)
chlam_norm=chlamydia_tot*100/(chlamydia_tot+proteobact_tot+mc_sizes)
proteobact_norm=proteobact_tot*100/(proteobact_tot+chlamydia_tot+mc_sizes)

##Figure 4A - metacell bacterial signal
pdf("bacterial_signal_per_mc.pdf",h=8,w=8,useDingbats=F)
plot(proteobact_norm,chlam_norm,pch=20,cex=4,col=as.character(xboc_ct_info$mc_color))
text(proteobact_norm,chlam_norm,labels=names(proteobact_norm),col=ifelse(pmax(proteobact_norm,chlam_norm)>4,"black",alpha("black",0)))
dev.off()

##Figure 4B/C - sc bacterial signal 2D projections
mc2d=scdb_mc2d("2dproj_v3")
scp_plot_gene_2d_metacell_bacteria(mc2d,mc,chlam_norm,out_fn="Chlamydia_2d.png",plot_mc=F,log=F)
scp_plot_gene_2d_metacell_bacteria(mc2d,mc,out_fn="Proteobacteria_2d.png",proteobact_norm,plot_mc=F,log=F)

##Figure 4E - fraction of infected/bacterial-containing cells per animal
cell_animal=as.vector(mat@cell_metadata[names(mc@mc),"dataset"])
names(cell_animal)=names(mc@mc)
frac_proteobact_cells=table(cell_animal[proteobact_positive_cells])*100/table(cell_animal)
frac_chlam_cells=table(cell_animal[chlam_positive_cells])*100/table(cell_animal)

pdf("Frac_infected_cells_per_animal.pdf",h=10,w=6,useDingbats=F)
par(mfrow=c(2,1))
barplot(frac_proteobact_cells,ylab="% infected cells",col="gray30",ylim=c(0,10),main="Proteobacteria")
barplot(frac_chlam_cells,ylab="% infected cells",col="gray30",ylim=c(0,10),main="Chalmydia")
dev.off()
```

