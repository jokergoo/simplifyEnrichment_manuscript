setwd("~/manuscript/simplifyEnrichment/figures")

library(cola)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)

data(golub_cola)

set.seed(123)

res = golub_cola["ATC:skmeans"]
df = get_signatures(res, k = 3, row_km = 3, plot = FALSE)
cl = get_classes(res, k = 3)[, 1]
cl[cl == 1] = 4; cl[cl == 2] = 1; cl[cl == 4] = 2
m = get_matrix(res)[df$which_row, ]
m = t(scale(t(m)))

ht = Heatmap(m, name = "Scaled\nexpression",
	col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
	top_annotation = HeatmapAnnotation(
		cola_class = cl,
		"ALL/AML" = get_anno(res)[, 1],
		col = list(cola_class = cola:::brewer_pal_set2_col[1:3],
			       "ALL/AML" = get_anno_col(res)[[1]])),
	row_split = paste0("km", df$km), show_row_dend = FALSE, show_row_names = FALSE,
	column_split = cl, show_column_dend = FALSE, show_column_names = FALSE,
	column_title = qq("A) @{nrow(df)} signature genes under FDR < 0.05")
)
p1 = grid.grabExpr(draw(ht, merge_legend = TRUE, padding = unit(c(12, 2, 2, 2), "mm")))

library(hu6800.db)
x = hu6800ENTREZID
mapped_probes = mappedkeys(x)
id_mapping = unlist(as.list(x[mapped_probes]))

lt = functional_enrichment(res, k = 3, id_mapping = id_mapping)

ago = c(rownames(lt[[1]]), rownames(lt[[2]]), rownames(lt[[3]]))
ago = unique(ago)
pm = matrix(1, nrow = length(ago), ncol = 3)
rownames(pm) = ago
colnames(pm) = c("km1", "km2", "km3")
	              
pm[lt[[1]]$ID, 1] = lt[[1]]$p.adjust
pm[lt[[2]]$ID, 2] = lt[[2]]$p.adjust
pm[lt[[3]]$ID, 3] = lt[[3]]$p.adjust

fdr_cutoff = 0.01
pm = pm[apply(pm, 1, function(x) any(x < fdr_cutoff)), ]
all_go_id = rownames(pm)

library(simplifyEnrichment)
sim_mat = GO_similarity(all_go_id)
col_fun_p = colorRamp2(c(0, -log10(fdr_cutoff), 4), c("#A1D76A", "#F7F7F7", "#C51B7D"))
ht_fdr = Heatmap(-log10(pm), col = col_fun_p, name = "FDR",
	show_row_names = FALSE, cluster_columns = FALSE,
	border = "black",
	heatmap_legend_param = list(at = c(0, -log10(fdr_cutoff), 4), 
		labels = c("1", fdr_cutoff, "< 0.0001")),
	width = unit(1.5, "cm"), use_raster = TRUE)
p2 = grid.grabExpr(
	simplifyGO(sim_mat, ht_list = ht_fdr, word_cloud_grob_param = list(max_width = 80), 
		verbose = FALSE, min_term = round(nrow(sim_mat)*0.01), control = list(partition_fun = partition_by_kmeanspp),
		column_title = qq("B) @{nrow(sim_mat)} GO terms clustered by 'binary cut'")
	), width = 14*2/3, height = 6.5)

library(cowplot)

pdf("figure7.pdf", width = 14, height = 6.5)
print(plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 2)))
dev.off()

