suppressMessages(source("../scripts/helper.R"))
library("zoo")

# list of species to compare
com_list = list( c("Xboc","Isopu"), c("Xboc","Spur"), c("Xboc","Hmia"), c("Xboc","Nvec") )

order_xboc = c("Xboc_GUT","Xboc_EPIDERMIS","Xboc_CILIATED_CELLS","Xboc_PHARYNGEAL","Xboc_PIGMENT","Xboc_GLAND1","Xboc_GLAND2","Xboc_GLAND3","Xboc_GLAND4","Xboc_MUSCLE","Xboc_MUSCLE.LIKE","Xboc_SECRETORY_EPITHELIIUM_1","Xboc_SECRETORY_EPITHELIIUM_2","Xboc_SECRETORY_EPITHELIIUM_3","Xboc_SECRETORY_EPITHELIIUM_4","Xboc_unk1","Xboc_unk2","Xboc_unk3","Xboc_EXCRETORY","Xboc_unk4","Xboc_NEURON","Xboc_unk5","Xboc_unk6","Xboc_unk7","Xboc_unk8","Xboc_unk9")

scp_select_top_markers = function(matrix, matrix_thr = 1.5, n_top_markers = 20, n_markers_rollmean = 2) {
	
	markers = unique(as.vector(unlist( apply(matrix, 2, function(c) { names( head(sort(-c[ c >= matrix_thr ]), n = n_top_markers ) ) } ))))
	markers_order = order(apply(matrix[markers,], 1, function(r) which.max(rollmean(r, n_markers_rollmean) )))
	markers_ordered = markers [ markers_order ]
	return(markers_ordered)
	
}



for (n in 1:length(com_list)) {

	# log	
	spi = com_list[[n]][1]
	spj = com_list[[n]][2]
	message(sprintf("load %s v %s", spi, spj))
	
	# load
	com = read.table(sprintf("results_alignment_samap/samap.%s-%s.scores.cts.top100q.csv", spi, spj), sep = "\t", header = TRUE, row.names = 1)
	rownames(com) = colnames(com)
	com = com [ grepl(sprintf("^%s_", spi), rownames(com)), grepl(sprintf("^%s_", spj), colnames(com)) ]
	com = com [ !grepl("Undefined", rownames(com)), !grepl("Undefined", colnames(com)) ]
	com = com [ !grepl("Unknown", rownames(com)), !grepl("Unknown", colnames(com)) ]
	com = com [ !grepl("uncharacterized", rownames(com)), !grepl("uncharacterized", colnames(com)) ]
	com = com [ !grepl("contamination", rownames(com)), !grepl("contamination", colnames(com)) ]
	com = com [ !grepl("Nvec_nan", rownames(com)), !grepl("Nvec_nan", colnames(com)) ]
	com = com [ order_xboc , ]
	
	order_other = scp_select_top_markers(t(com), 0, 20, 1)
	com = com [ , order_other ]

	# print	
	message(sprintf("plot %s v %s", spi, spj))
	hm = plot_complex_heatmap(
		com,
		name = "samap",
		color_min = 0,
		color_max = 1,
		color_mat = c("white","#d6e72e","#6fb600","#003f4d"),
		cluster_row = FALSE,
		cluster_col = FALSE,
		do_dotplot = TRUE,
		use_raster = FALSE,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		cex_dotplot = 0.03
	)
	
	pdf(sprintf("results_alignment_samap/samap.%s-%s.scores.cts.top100q.pdf", spi, spj), height = 6, width = 6)
	print(hm)
	dev.off()
	
}