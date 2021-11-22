setwd("~/manuscript/simplifyEnrichment/figures")

library(ComplexHeatmap)
library(cowplot)
library(circlize)
library(GetoptLong)
library(ggplot2)
library(simplifyEnrichment)

load("compare_similarity_random_BP.RData")

df1 = do.call(rbind, lapply(n_set_list, function(x) data.frame(x = x, cate = names(x))))
df1$type = "All sizes"
df2 = do.call(rbind, lapply(large_n_set_list, function(x) data.frame(x = x, cate = names(x))))
df2$type = "Size >= 5"
df3 = do.call(rbind, lapply(avg_n_set_list, function(x) data.frame(x = x, cate = names(x))))
df3$type = "All sizes"
df4 = do.call(rbind, lapply(large_avg_n_set_list, function(x) data.frame(x = x, cate = names(x))))
df4$type = "Size >= 5"

p1 = ggplot(rbind(df1, df2), 
	aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x, col = type)) + 
	geom_boxplot() + labs(x = "", y = "Numbers of clusters") + 
	scale_y_log10(limits = c(1, 500), breaks = c(1, 10, 100, 500)) + 
	theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
	ggtitle("A) On random GO lists")
p2 = ggplot(rbind(df3, df4), 
	aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x, col = type)) + 
	geom_boxplot() + labs(x = "", y = "Average numbers of terms per cluster") + 
	scale_y_log10(limits = c(1, 500), breaks = c(1, 10, 100, 500)) + 
	theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
	ggtitle("B) On random GO lists")



arr = array(dim = c(dim(con_mat_list[[1]]), length(con_mat_list)))
for(i in seq_along(con_mat_list)) {
	arr[, , i] = con_mat_list[[i]]
}
con_mat = apply(arr, 1:2, mean)
dimnames(con_mat) = dimnames(con_mat_list[[1]])

con_col_fun = colorRamp2(c(0, 1), c("white", "orange"))
ht_opt$TITLE_PADDING = unit(2.5, "mm")
p3 = grid.grabExpr(draw(Heatmap(con_mat, show_heatmap_legend = FALSE,
	col = con_col_fun, cell_fun = function(j, i, x, y, w, h, fill) {
		grid.text(sprintf("%.2f", con_mat[i, j]), x, y, gp = gpar(fontsize = 8))
	}, column_title = "C) Mean concordance",
	row_dend_width = unit(5, "mm"), column_dend_height = unit(5, "mm"),
	row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8)),
	padding = unit(c(4, 4, 0, 1), "mm")))


clt1 = readRDS("clt_random_BP_10.rds")
p4 = grid.grabExpr(ht_clusters(clt1[[1]], clt1[[2]]$binary_cut,  draw_word_cloud = FALSE, 
	col = colorRamp2(c(0, quantile(clt1[[1]], 0.95)), c("white", "red")),
	column_title = "D) semantic similarity", show_heatmap_legend = FALSE, padding = unit(c(10, 4, 0, 1), "mm")))
clt2 = readRDS("clt_random_BP_kappa_10.rds")
p5 = grid.grabExpr(ht_clusters(clt2[[1]], clt2[[2]]$binary_cut,  draw_word_cloud = FALSE, 
	col = colorRamp2(c(0, quantile(clt2[[1]], 0.95)), c("white", "red")),
	column_title = "E) kappa similarity", show_heatmap_legend = FALSE, padding = unit(c(10, 4, 0, 1), "mm")))



lgd_list = list(
	Legend(title = "Numbers of clusters", labels = c("All sizes", "Size >= 5"), type = "boxplot", legend_gp = gpar(col = c("#F8766D", "#00BFC4"), fill = "white"), direction = "horizontal", nrow = 1),
	Legend(title = "Mean concordance", at = c(0, 0.25, 0.5, 0.75, 1), col_fun = con_col_fun, direction = "horizontal"),
	Legend(title = " ", at = "All other small\nclusters (size < 5)", legend_gp = gpar(fill = "darkgreen"), grid_height = unit(8, "mm"), grid_width = unit(2, "mm"), direction = "horizontal"),
	Legend(title = "Semantic similarity\nrandom GO lists", col_fun = colorRamp2(c(0, quantile(clt1[[1]], 0.95)), c("white", "red")), direction = "horizontal", legend_width = unit(4.5, "cm")),
	Legend(title = "Kappa similarity\nrandom GO lists", col_fun = colorRamp2(c(0, quantile(clt2[[1]], 0.95)), c("white", "red")), at = c(0, 0.005, 0.01, 0.015, 0.02), direction = "horizontal", legend_width = unit(4.5, "cm"))
)
###########

load("compare_similarity_EBI_Expression_Atlas_GO_BP.RData")

df1 = do.call(rbind, lapply(n_set_list, function(x) data.frame(x = x, cate = names(x))))
df1$type = "All sizes"
df2 = do.call(rbind, lapply(large_n_set_list, function(x) data.frame(x = x, cate = names(x))))
df2$type = "Size >= 5"
df3 = do.call(rbind, lapply(avg_n_set_list, function(x) data.frame(x = x, cate = names(x))))
df3$type = "All sizes"
df4 = do.call(rbind, lapply(large_avg_n_set_list, function(x) data.frame(x = x, cate = names(x))))
df4$type = "Size >= 5"

p6 = ggplot(rbind(df1, df2), 
	aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x, col = type)) + 
	geom_boxplot() + labs(x = "", y = "Numbers of clusters") + 
	scale_y_log10(limits = c(1, 500), breaks = c(1, 10, 100, 500)) + 
	theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
	ggtitle("F) On Expression Atlas datasets")
p7 = ggplot(rbind(df3, df4), 
	aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x, col = type)) + 
	geom_boxplot() + labs(x = "", y = "Average numbers of terms per cluster") + 
	scale_y_log10(limits = c(1, 500), breaks = c(1, 10, 100, 500)) + 
	theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
	ggtitle("G) On Expression Atlas datasets")


arr = array(dim = c(dim(con_mat_list[[1]]), length(con_mat_list)))
for(i in seq_along(con_mat_list)) {
	arr[, , i] = con_mat_list[[i]]
}
con_mat = apply(arr, 1:2, mean)
dimnames(con_mat) = dimnames(con_mat_list[[1]])

con_col_fun = colorRamp2(c(0, 1), c("white", "orange"))
ht_opt$TITLE_PADDING = unit(2.5, "mm")
p8 = grid.grabExpr(draw(Heatmap(con_mat, show_heatmap_legend = FALSE,
	col = con_col_fun, cell_fun = function(j, i, x, y, w, h, fill) {
		grid.text(sprintf("%.2f", con_mat[i, j]), x, y, gp = gpar(fontsize = 8))
	}, column_title = "H) On Expression Atlas datasets",
	row_dend_width = unit(5, "mm"), column_dend_height = unit(5, "mm"),
	row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8)),
	padding = unit(c(4, 4, 0, 1), "mm")))


ht_opt$TITLE_PADDING = unit(2.5, "mm")
clt1 = readRDS("clt_E-GEOD-10718_A-AFFY-44_g4_g3_GO_BP.rds")
p9 = grid.grabExpr(ht_clusters(clt1[[1]], clt1[[2]]$binary_cut,  draw_word_cloud = FALSE, 
	col = colorRamp2(c(0, quantile(clt1[[1]], 0.95)), c("white", "red")),
	column_title = "I) semantic similarity", show_heatmap_legend = FALSE, padding = unit(c(10, 4, 0, 1), "mm")))
clt2 = readRDS("clt_E-GEOD-10718_A-AFFY-44_g4_g3_GO_BP_kappa.rds")
p10 = grid.grabExpr(ht_clusters(clt2[[1]], clt2[[2]]$binary_cut,  draw_word_cloud = FALSE, 
	col = colorRamp2(c(0, quantile(clt2[[1]], 0.95)), c("white", "red")),
	column_title = "J) kappa similarity", show_heatmap_legend = FALSE, padding = unit(c(10, 4, 0, 1), "mm")))


lgd_list = c(
	lgd_list,
	list(
		Legend(title = "Semantic similarity\nExpression Atlas datasets", col_fun = colorRamp2(c(0, quantile(clt1[[1]], 0.95)), c("white", "red")), direction = "horizontal", legend_width = unit(4.5, "cm")),
		Legend(title = "Kappa similarity\nExpression Atlas datasets", col_fun = colorRamp2(c(0, quantile(clt2[[1]], 0.95)), c("white", "red")), direction = "horizontal", legend_width = unit(4.5, "cm"))
	)
)

df = do.call(rbind, lapply(n_set_list, function(x) data.frame(x = x, cate = names(x))))
df$cate = factor(df$cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa"))
p11 = ggplot(df, aes(x = x, fill = cate)) + geom_histogram() + facet_grid(~cate) + theme(legend.position = "none") +
	labs(x = "Numbers of clusters", y = "Counts") + 
	scale_x_continuous(breaks = c(1, 50, 100, 150, 200)) +
	ggtitle("K) Number of clusters, on Expression Atlas datasets")

df = do.call(rbind, lapply(max_cluster_prop_list, function(x) data.frame(x = x, cate = names(x))))
df$cate = factor(df$cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa"))
p12 = ggplot(df, aes(x = x, fill = cate)) + geom_histogram() + facet_grid(~cate) + theme(legend.position = "none") +
	labs(x = "Fraction of the largest cluster", y = "Counts") + ggtitle("L) Fraction of the largest cluster, on Expression Atlas datasets")


lgd = packLegend(list = lgd_list, direction = "horizontal", column_gap = unit(6, "mm"))
lgd_h = convertHeight(grobHeight(lgd@grob), "inch", valueOnly = TRUE)*1.1
pdf("figure6.pdf", width = 14, height = 14)
print(plot_grid(
	plot_grid(p1, p2, p3, p4, p5, 
	            p6, p7, p8, p9, p10, 
	            nrow = 2),
	lgd@grob, p11, p12, nrow = 4, rel_heights = c((14-lgd_h)/2, lgd_h, (14-lgd_h)/4, (14-lgd_h)/4))
)
dev.off()
