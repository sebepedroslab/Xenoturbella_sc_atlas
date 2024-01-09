###Tested with R version 3.5.
library(metacell)
library(plyr)
library(reshape2)
library(dplyr)
library(zoo)
library('RColorBrewer')
library(scales)
library(pheatmap)
library(ggplot2)


scr_load_cell_type_table=function(input_table,mc_object){
  
  cell_type_table=read.table(input_table,h=TRUE,sep="\t",comment.char="",stringsAsFactors=F)
  rownames(cell_type_table)=cell_type_table$metacell
  
  #niche_order<<-as.character(scr_cell_type_table$metacell)
  sc_ct_label=as.vector(cell_type_table[as.character(mc_object@mc),"cell_type"])
  names(sc_ct_label)=names(mc_object@mc)
  
  clust_cols=as.character(cell_type_table[,"color"])
  names(clust_cols)=rownames(cell_type_table)
  
  cells_cols=cell_type_table[as.character(mc_object@mc),"color"]
  cells_cols=as.character(cells_cols)
  names(cells_cols)=names(mc_object@mc)
  
  mc_object@colors=clust_cols
  return(list(ct_table=cell_type_table,mc_color=clust_cols,sc_ct_label=sc_ct_label,sc_color=cells_cols))
}


scr_cell_type_fp=function(input_table,mc_object,mat_object){  ##compute cell_tep footprint based on a standard metacell->cell type definition table, using same strategy as mc_fp.
	cell_type_table=read.table(input_table,h=TRUE,sep="\t",comment.char="")
	rownames(cell_type_table)=cell_type_table$metacell
	
	sc_ct_label=as.vector(cell_type_table[as.character(mc_object@mc),"cell_type"])
	names(sc_ct_label)=names(mc_object@mc)
	
	cells_cols=cell_type_table[as.character(mc_object@mc),"color"]
	cells_cols=as.character(cells_cols)
	names(cells_cols)=names(mc_object@mc)
	
	umis=as.matrix(mat_object@mat)
	umis=umis[rownames(mc_object@mc_fp),names(mc_object@mc)] #filter low expression genes not included in the mc_fp

	
	ct_geomean=t(apply(umis, 1,  function(x) tapply(x, sc_ct_label,function(y) exp(mean(log(1+y)))-1)))
	ct_meansize=tapply(colSums(umis), sc_ct_label, mean)
	ideal_cell_size=pmin(1000,median(ct_meansize))
	g_fp=t(ideal_cell_size * t(ct_geomean)/as.vector(ct_meansize))
	fp_reg=0.05
	g_fp_n=(fp_reg + g_fp)/apply(fp_reg + g_fp, 1, median)

	##for compatibility with other functions, return as a MC-like object
	ct_table=mc_object
	ct_table@mc_fp=g_fp_n[,as.character(unique(cell_type_table$cell_type))]
	ct_table@mc=sc_ct_label
	ct_table@colors=cells_cols
	ct_table@cell_names=names(sc_ct_label)
	return(ct_table)
}


scr_recompute_fp=function(cell_classification,mc_object,mat_object){ #cell_classification is a named vector classifying cells. Returns a MC object
  															#E.g. usefull to examine a Seurat or scanpy classification.
	
	cluster_cols=colorRampPalette(c("darkgray", "burlywood1", "chocolate4","orange", "red", "purple", "blue","darkgoldenrod3", "cyan"))(length(unique(cell_classification)))
	
	umis=as.matrix(mat_object@mat)
	umis=umis[,names(cell_classification)]
	ct_geomean=t(apply(umis, 1,  function(x) tapply(x, cell_classification,function(y) exp(mean(log(1+y)))-1)))
	ct_meansize=tapply(colSums(umis), cell_classification, mean)
	ideal_cell_size=pmin(1000,median(ct_meansize))
	g_fp=t(ideal_cell_size * t(ct_geomean)/as.vector(ct_meansize))
	fp_reg=0.05
	g_fp_n=(fp_reg + g_fp)/apply(fp_reg + g_fp, 1, median)
	##for compatibility with other functions, return as a MC-like object
	ct_table=mc_object
	ct_table@mc_fp=g_fp_n
	ct_table@mc=cell_classification
	ct_table@colors=cluster_cols
	ct_table@cell_names=names(cell_classification)
	return(ct_table)
}

scr_plot_cmod_markers_ct_colors = function(mc_object,mat_object,output_file,gene_annot_file,sn_table=NULL,mc_object_gene_selection=NULL,black_list=c(),genes_to_add=c(),
										height=60, width=60,transversality_N=ncol(mc_object@mc_fp),per_clust_genes=25,
										gene_min_fold=2,clust_ord=NULL,lateral_colorbar=FALSE,plot_sc=FALSE,clust_color=NULL,p_min=5)
  {    
    
	if(is.null(clust_ord)){ clust_ord=colnames(mc_object@mc_fp) }
	if(is.null(sn_table)){lateral_colorbar=FALSE} #double-check, there can be no possible gene-CT assignment if no cell type table is provided
	if(is.null(sn_table) & is.null(clust_color)){
		stop("...You must provide either a color vector or a cell type annotation table...")
	}
	
	if(!is.null(sn_table)){
		cell_type_table=read.table(sn_table,h=TRUE,sep="\t",comment.char="")
		rownames(cell_type_table)=cell_type_table$metacell
	}
	
	annot=read.table(gene_annot_file,header=T,sep="\t",fill=TRUE,quote="",row.names=1)
  
	niche_geomean_n=as.matrix(mc_object@mc_fp[!grepl("peak",rownames(mc_object@mc_fp)),])
	
	if(is.null(mc_object_gene_selection)){
		genes=unique(as.vector(unlist(apply(niche_geomean_n, 2, function(x) names(head(sort(-x[x>gene_min_fold]),n=per_clust_genes))))))
		transversal_genes=names(which(apply(niche_geomean_n, 1, function(x) sort(x,decreasing=T)[transversality_N]>1.8)))
		genes=setdiff(genes, transversal_genes)
	}else{
		genes=unique(as.vector(unlist(apply(mc_object_gene_selection@mc_fp, 2, function(x) names(head(sort(-x[x>gene_min_fold]),n=per_clust_genes))))))
	}
	
	genes=setdiff(genes, black_list)
	genes=union(genes,genes_to_add)
	genes=intersect(genes,rownames(niche_geomean_n))
	
	message(length(genes))
	
	
	message("Will use ",length(genes)," genes")
	
	mat_niche=as.matrix(niche_geomean_n[genes,])
    #gene_ord=order(apply(mat_niche[,clust_ord],1,function(x) which.max(rollmean(x,1))))
	genes_reord=genes[order(apply(mat_niche[,clust_ord],1,function(x) which.max(rollmean(x,1))))]
	

	if(is.null(clust_color)){

		clust_color=as.character(cell_type_table[as.character(clust_ord),"color"])
	}
	
    pdf(paste0(output_file,"_mc.pdf"), h=height, w=width,useDingbats=F)
    par(mar=c(0,0,0,0))
    par(fig=c(0.2,0.8,0.1,0.9))
    shades=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
    #image(t(pmax(log2(mat_niche_to_plot[, clust_ord]),0)), col=shades,xaxt="n",yaxt="n")
    image(t(pmin(log2(niche_geomean_n[genes_reord,as.character(clust_ord)]+1),p_min)), col=shades,xaxt="n",yaxt="n")
	
	if(lateral_colorbar==FALSE){
		mtext(annot[genes_reord,2], side=4, at=seq(0,1,length.out=length(genes_reord)), col=ifelse(genes_reord %in% genes_to_add,"darkred","black"),las=1, line=1)
    }
          
    mtext(paste(annot[genes_reord,1],genes_reord,sep="||"), side=2,col=ifelse(genes_reord %in% genes_to_add,"darkred","black"), at=seq(0,1,length.out=length(genes_reord)), las=1, line=1)
    mtext(clust_ord,side=1,at=seq(0,1,length.out=sum(ncol(mat_niche))), las=2,line=3,cex=3)

    mtext(clust_ord,side=3,at=seq(0,1,length.out=sum(ncol(mat_niche))), las=2,line=3,cex=3)

	par(fig=c(0.2,0.8,0.03,0.06),new=TRUE)
	image(as.matrix(1:length(clust_ord)),col=clust_color, axes = F,xaxt='n',yaxt='n')
	
	if(lateral_colorbar==TRUE){
		gene_to_metaclust=as.character(clust_ord[apply(mat_niche[genes_reord,clust_ord],1,function(x) which.max(rollmean(x,1)))])
		par(fig=c(0.82,0.85,0.1,0.9),new=TRUE)
		image(t(as.matrix(1:length(gene_to_metaclust))),col=as.character(cell_type_table[gene_to_metaclust,"color"]), axes = F,xaxt='n',yaxt='n')
	}
    
    dev.off()
	
	print(length(genes))
	
	if(plot_sc==TRUE){
	##################PLOT SINGLE-CELL PROFILE########################
	cell_order=c()  
    for (niche in as.character(clust_ord)){
      cells=names(mc_object@mc[which(mc_object@mc==niche)])    
	  cell_order=c(cell_order,cells)
    }
	cluster_cell_count=as.matrix(table(mc_object@mc))
    n_cells_cluster=cluster_cell_count[clust_ord,1]

	mat = as.matrix(mat_object@mat[genes, cell_order])
	totu = colSums(as.matrix(mat_object@mat[, cell_order]))
	mat = t(t(mat)/totu)*800

	lus_1 = log2(1+7*mat[genes[gene_ord], cell_order])
	lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
	lus_smoo = t(apply(lus[genes[gene_ord],cell_order], 1, function(x) rollmean(x,4, fill=0)))

    pdf(paste0(output_file,"_sc_cells.pdf"), h=height, w=width*2,useDingbats=F)
    par(mar=c(0,0,0,0))
    par(fig=c(0.2,0.90,0.1,0.9))
    shades=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
    image(t(pmin(lus_smoo,3.8)), col=shades,xaxt="n",yaxt="n")
	#print(quantile(lus_smoo,seq(0,1,by=0.01)))
    x=0
    for(i in 1:length(n_cells_cluster)) {
      abline(v=(n_cells_cluster[i]+x)/(sum(n_cells_cluster)-1)-1/(2*sum(n_cells_cluster)), lwd=4,col="grey")      
      mtext(clust_ord[i], side=1, at=((n_cells_cluster[i]/1.5+x)/(sum(n_cells_cluster)-1)), adj=1, las=2, line=1,cex=5)
      mtext(clust_ord[i], side=3, at=((n_cells_cluster[i]/1.5+x)/(sum(n_cells_cluster)-1)), las=2, line=1,cex=5)     
      x=x+n_cells_cluster[i]
    }
	par(fig=c(0.2,0.90,0.02,0.05),new=TRUE)

	image(as.matrix(1:length(cell_order)),
        col=as.character(cell_type_table[as.character(mc_object@mc[cell_order]),"color"]), axes = F,xaxt='n',yaxt='n')
	
	gene_to_metaclust=as.character(clust_ord[apply(mat_niche[gene_ord,clust_ord],1,function(x) which.max(rollmean(x,1)))])
	par(fig=c(0.92,0.95,0.1,0.9),new=TRUE)
	image(t(as.matrix(1:length(gene_to_metaclust))),col=as.character(cell_type_table[gene_to_metaclust,"color"]), axes = F,xaxt='n',yaxt='n')
	
	x=sapply(clust_ord,function(x) length(which(mc_object@mc_fp[,x] > 2 & apply(mc_object@mc_fp[,setdiff(clust_ord,x)],1,max)<2)))
	barwidths=table(gene_to_metaclust)[clust_ord]
	barwidths[is.na(barwidths)]<-0
	par(fig=c(0.11,0.19,0.1,0.9),new=TRUE)

	barplot(-log2(x+1),barwidths,col=as.character(cell_type_table$color),cex.axis=4,xaxt='n',yaxt='n',horiz=TRUE,xaxs='i',yaxs='i',space=0.05)
	axis(1,at=-0:-max(log2(x+1)),labels=0:max(log2(x+1)),cex=6,lwd=2,line=1,cex.lab=3)
	dev.off()
	}
}


scr_dot_plot_map=function(umi_frac_object,mc_object,markers_file,ct_info,out_fn,min_gene_fc=1,sort_genes=F,width=20,height=25,fc_ths=3,ct_mode=F){
	markers_table=read.table(markers_file,header=FALSE,row.names=1,sep="\t")
	markers=intersect(rownames(markers_table),rownames(mc_object@mc_fp))
	f_marker_fc=apply(mc_object@mc_fp[markers,],1,max) > min_gene_fc
	markers=intersect(markers,names(which(f_marker_fc)))
	
	if(sort_genes==TRUE){markers_sorted=markers[as.numeric(order(apply(mc_object@mc_fp[markers,], 1,function(x) which.max(rollmean(x,1)))))]} else {markers_sorted=rev(markers)}
	
	mat_fc=pmin(mc_object@mc_fp[markers_sorted,],fc_ths)
	mat_umifrac=as.matrix(umi_frac_object[markers_sorted,])
	mat_umifrac=pmin(mat_umifrac,quantile(mat_umifrac,0.98))
	base_mat=matrix(ncol=ncol(mat_umifrac),nrow=nrow(mat_umifrac))
	base_mat[is.na(base_mat)]=1
	colnames(base_mat)=colnames(mat_umifrac) #####<---
	
	base_mat=melt(base_mat)
	colnames(base_mat)=c("genes","metacells","na")
	base_mat$umifrac=melt(mat_umifrac)$value
	base_mat$fc=melt(mat_fc)$value
	if(ct_mode==FALSE){
		base_mat$ct=ct_info$ct_table[base_mat$metacells,"cell_type"]
	}else{
		base_mat$ct=base_mat$metacells
	}
	
	ct_color=as.character(unique(ct_info$ct_table$color))
	names(ct_color)=as.character(unique(ct_info$ct_table$cell_type))

	b=ggplot(base_mat, aes(x=metacells,y=genes)) + geom_point(aes(size=umifrac,colour=ct)) + scale_colour_manual(values=ct_color) + 
			scale_y_continuous(breaks=1:length(markers_sorted),labels=markers_sorted, sec.axis=dup_axis(labels=substring(as.character(markers_table[markers_sorted,1]),1,70))) + theme(axis.text.y = element_text(angle=0)) +
			theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
					
	a=ggplot(base_mat, aes(x=metacells,y=genes)) + geom_point(aes(size=umifrac,colour=fc)) + scale_colour_gradientn(colors=c("lightgrey","lightgrey","orange","red","purple","black")) + 
			scale_y_continuous(breaks=1:length(markers_sorted),labels=markers_sorted, sec.axis=dup_axis(labels=substr(as.character(markers_table[markers_sorted,1]),1,70))) + theme(axis.text.y = element_text(angle=0)) +
			theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
				
	d=ggplot(base_mat, aes(x=metacells,y=genes)) + geom_point(aes(size=umifrac,alpha=fc,colour=ct)) + scale_colour_manual(values=ct_color) + 
			scale_y_continuous(breaks=1:length(markers_sorted),labels=markers_sorted, sec.axis=dup_axis(labels=substring(as.character(markers_table[markers_sorted,1]),1,70)))  + theme(axis.text.y = element_text(angle=0)) +
			theme(axis.line = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
			
	
	#pdf(file="Dotmap_fc_umifrac.pdf",h=height,w=width,useDingbats=F);print(a);dev.off()
	#pdf(file="Dotmapumifrac.pdf",h=height,w=width,useDingbats=F);print(b);dev.off()
	pdf(file=paste0(out_fn,"_dotmapumifrac_fc_alpha.pdf"),h=height,w=width,useDingbats=F);print(d);dev.off()

}


scr_barplot_heatmap_markers=function(mc_object,mat_object,mc_counts,markers_file,barplot_dir=".",heatmap_file,w=30,
                                     T_totumi=15,clust_order=NULL,
                                     mc_color=NULL,
                                     print_barplots=FALSE,pmin=3,
                                     min_gene_fc=2,set_pmin=FALSE,column_names=NULL,write_table=FALSE,sort_markers=TRUE){
  if(is.null(clust_order)){
  	clust_order=colnames(mc_object@mc_fp)
  }
  
  clusts=mc_object@mc
  footprint_table=mc_object@mc_fp[,clust_order]

  if(!is.null(mc_color)){
    color=as.character(mc_color)
  }else { color="gray33"}
	
  markers_table=read.table(markers_file,header=FALSE,row.names=1,sep="\t")
  markers=intersect(rownames(markers_table),rownames(footprint_table))
  
  f_marker_cov = rowSums(mc_counts[markers,]) > T_totumi
  f_marker_fc=apply(mc_object@mc_fp[markers,],1,max) > min_gene_fc
  
  
  
  markers=intersect(names(which(f_marker_cov)),names(which(f_marker_fc)))
  markers=intersect(markers,rownames(footprint_table))
  ###define heatmap height based on number of markers
  h=pmin(pmax(length(markers)*0.7,10),100)
  #h=pmin(pmax(length(markers)*60,1500),12000)
  
  
  #marker_fp = scr_niche_footprint[markers,as.character(niche_order)]
  #marker_fp = scr_niche_umifrac[markers,as.character(niche_order)]
  marker_fp = footprint_table[markers,as.character(clust_order)]
  if(sort_markers){ markers_sorted=markers[as.numeric(order(apply(marker_fp, 1,function(x) which.max(rollmean(x,1)))))]} else {markers_sorted=rev(markers)}

  
  markers_tot_umi=rowSums(mc_counts[markers_sorted,])
  marker_max_umi_frac=round(apply(mc_counts[markers_sorted,]/markers_tot_umi,1,max),2)
  
  
  print(paste0(length(markers_sorted),"  markers survived"))
  shades2=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)

  if(set_pmin){
	pmin=pmin(quantile(footprint_table[markers_sorted,clust_order],0.995),4)
  }
  message("I'm using PMIN   ", pmin)
  marker_fp_to_plot=pmin(footprint_table[markers_sorted,clust_order],pmin)
  
  #marker_fp_to_plot=pmax(pmin(log2(footprint_table[markers_sorted,niche_order]),2),-2)
	pdf(heatmap_file,height=h,width=w,useDingbats=F) 
	par(mar=c(0,0,0,0))
	par(fig=c(0.3,0.7,0.2,0.95))
    image(t(marker_fp_to_plot[markers_sorted,]), xaxt='n',yaxt='n', col=shades2)
    if(is.null(column_names)){
		if(is.null(mc_color)){mtext(clust_order, side=1,at=seq(0,1,length.out=length(clust_order)), las=1,cex=2,line=0.2,las=2)}
		mtext(clust_order, side=3,at=seq(0,1,length.out=length(clust_order)), las=1,cex=1,line=0.2,las=2)
	} else {
		if(is.null(mc_color)){mtext(column_names, side=1,at=seq(0,1,length.out=length(column_names)), las=1,cex=2,line=0.2,las=2)}
		mtext(column_names, side=3,at=seq(0,1,length.out=length(column_names)), las=1,cex=1,line=0.2,las=2)
	}
	if(ncol(markers_table)==1){
		mtext(paste(markers_sorted,markers_tot_umi,marker_max_umi_frac,sep=" "), line=2,side=2,at=seq(0,1,length.out=length(markers_sorted)), las=1,cex=2,adj=1)
	}else{
		mtext(paste(markers_table[markers_sorted,2],markers_sorted,sep="||"), side=2,at=seq(0,1,length.out=length(markers_sorted)), line=2,las=1,cex=2,adj=1)
	}
    
    mtext(markers_table[markers_sorted,1], side=4,at=seq(0,1,length.out=length(markers_sorted)), line=2,las=1,cex=2,adj=0)
    

	for(i in 1:length(clust_order)){       abline(v=(i-0.5)/(length(clust_order)-1), lwd=0.5)  }
	for(i in 1:length(markers_sorted)) {       abline(h=(i-0.5)/(length(markers_sorted)-1), lwd=0.5)   }  
      
	if(!is.null(mc_color)){
		par(fig=c(0.3,0.70,0.15,0.18),new=TRUE)
		image(as.matrix(1:length(mc_color)), col=as.character(mc_color), axes = F,xaxt='n',yaxt='n')
	}
	dev.off()  
	
	if(write_table){
		markers_table=as.matrix(markers_table)
		m_to_write=cbind(round(footprint_table[rev(markers_sorted),clust_order],2),as.character(markers_table[rev(markers_sorted),1]),as.character(markers_table[rev(markers_sorted),2]))
		colnames(m_to_write)=c(clust_order,"BBH","Pfam_domains")
		write.table(m_to_write,file=paste0(heatmap_file,"_table.txt"),col.names=T,row.names=T,sep="\t",quote=F)
	}
  
	if(print_barplots){
		for(nm in markers) {
		fn = sprintf("%s/%s_barplot.png", barplot_dir, nm)
		png(fn,h=300,w=1500)
		par(mar=c(3,8,3,3))
		#      barplot(scr_niche_umifrac[nm,],main=markers_table[nm,1],col=as.character(scr_cell_type_table[niche_order,"color"]),cex.axis=3,cex.main=2,las=2)
		barplot(scr_niche_umifrac[nm,niche_order],main=markers_table[nm,1],cex.axis=3,cex.main=2,las=2,border=F,col=color,space=0.1)
		dev.off()
		} 
	}
}


scp_plot_gene_2d_metacell_bacteria = function(sc2d_object,mc_object,norm_bacteria_exp,out_fn="Bact",plot_mc=F,log=T,quant_pmin=0.95){ 
	
	#define color per mc
	if(log==T){expr=log2(norm_bacteria_exp)} else { expr= pmin(norm_bacteria_exp,quantile(norm_bacteria_exp,quant_pmin))}
	
	
	cols=colorRampPalette(brewer.pal(7,"PuBuGn"))(1000)[as.numeric(cut(expr,breaks=1000))]
	#cols=viridis(1000)[as.numeric(cut(expr,breaks=1000))]
	#cols=colorRampPalette(c("white","grey90","grey50","black","brown4"))(1000)[as.numeric(cut(expr,breaks=1000))]

	#cols=colorRampPalette(c("lightgray", "orange","red","purple","black"))(100)[as.numeric(cut(expr,breaks=100))]

	png(paste0(out_fn,".png"),h=2500,w=2500)
	par(fig=c(0.1,0.95,0.1,0.95))
	cell_colors=rep(cols,table(mc_object@mc))
	plot(sc2d_object@sc_x[names(sort(mc_object@mc))],sc2d_object@sc_y[names(sort(mc_object@mc))],col=alpha(cell_colors,0.4),pch=20,cex=5,yaxt="n",xaxt="n",xlab="",ylab="")
	if(plot_mc){
	points(sc2d_object@mc_x,sc2d_object@mc_y,col="black",pch=21,cex=10,bg=cols,lwd=1)
	text(sc2d_object@mc_x,sc2d_object@mc_y,cex=3)
	}
	par(fig=c(0.1,0.3,0.03,0.11),new=TRUE)   ####FC legend
		x=seq(min(expr),max(expr),(max(expr))/(length(cols)-1))
	col_legend=colorRampPalette(brewer.pal(7,"PuBuGn"))(1000)
	#col_legend=colorRampPalette(c("grey90","grey90","grey50","brown4"))(1000)
	
	image(x=x, y=c(0,1), col=col_legend, yaxt="n", z=matrix(nrow=length(x),ncol=1,data=c(1:length(x))),
        ylab="",xlab="",cex.axis=3,xaxt="n")
	mtext(side=2,round(min(expr),1),las=2,cex=3)
	mtext(side=4,round(max(expr),1),las=2,cex=3)
	mtext(side=1,"Bacterial signal",cex=4,line=2)
	dev.off()

}


scp_plot_gene_2d_sc_bacteria = function(sc2d_object,mc_object,bacteria_exp,log=T,quant_pmin=0.95,out_fn="Bact"){ 
	
	#define color per mc
	if(log==T){expr=log2(bacteria_exp+1)} else { expr= pmin(bacteria_exp,quantile(bacteria_exp,quant_pmin))}
	
	cols=colorRampPalette(brewer.pal(7,"PuBuGn"))(1000)[as.numeric(cut(expr,breaks=1000))]
	#cols=viridis(1000)[as.numeric(cut(expr,breaks=1000))]
	cols=colorRampPalette(c("grey80","grey70","grey40","black","brown4"))(1000)[as.numeric(cut(expr,breaks=1000))]
	names(cols)=names(expr)
	#cols=colorRampPalette(c("lightgray", "orange","red","purple","black"))(100)[as.numeric(cut(expr,breaks=100))]

	png(paste0(out_fn,".png"),h=2500,w=2500)
	par(fig=c(0.1,0.95,0.1,0.95))
	cell_colors=cols[names(sort(mc_object@mc))]
	plot(sc2d_object@sc_x[names(sort(mc_object@mc))],sc2d_object@sc_y[names(sort(mc_object@mc))],col=alpha(cell_colors,0.5),pch=20,cex=c(2,5,8,11)[as.numeric(cut(expr[names(sort(mc_object@mc))],breaks=4))],yaxt="n",xaxt="n",xlab="",ylab="")

	par(fig=c(0.1,0.3,0.03,0.11),new=TRUE)   ####FC legend
		x=seq(min(expr),max(expr),(max(expr))/(length(cols)-1))
	col_legend=colorRampPalette(c("white","grey90","grey50","black","brown4"))(1000)
	#col_legend=colorRampPalette(c("grey90","grey90","grey50","brown4"))(1000)

	image(x=x, y=c(0,1), col=col_legend, yaxt="n", z=matrix(nrow=length(x),ncol=1,data=c(1:length(x))),
        ylab="",xlab="",cex.axis=3,xaxt="n")
	mtext(side=2,round(min(expr),1),las=2,cex=3)
	mtext(side=4,round(max(expr),1),las=2,cex=3)
	mtext(side=1,"Bacterial signal",cex=4,line=2)
	dev.off()

}
