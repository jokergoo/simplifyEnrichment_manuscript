library(simplifyEnrichment)

setwd("~/manuscript/simplifyEnrichment/figures")

set.seed(1234)
go_id = random_GO(500)
mat = GO_similarity(go_id)
cl = binary_cut(mat)

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.5), c("white", "red"))

cl = binary_cut(mat)
tb = sort(table(cl), decreasing = TRUE)
ind1 = which(as.character(cl) %in% names(tb)[1])
ind2 = which(as.character(cl) %in% names(tb)[1:2])

mat1 = mat[ind1, ind1]
km = kmeans(mat1, centers = 2)$cluster
plot1 = grid.grabExpr({
	draw(Heatmap(mat1, name = "Similarity", col = col_fun, show_row_dend = FALSE, show_column_dend = FALSE, 
		row_split = km, column_split = km,
		show_row_names = FALSE, show_column_names = FALSE,
		border = "#404040", row_title = NULL, 
		column_title = gt_render("A) Similarity matrix <span style='font-family:Times;'>**M**<sub>*a*</sub></span>"),
		row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_heatmap_legend = TRUE))
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>11</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 1)
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>12</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 2)
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>21</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 1)
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>22</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 2)
})
mat2 = mat[ind2, ind2]
km = kmeans(mat2, centers = 2)$cluster
plot2 = grid.grabExpr({
	draw(Heatmap(mat2, name = "Similarity", col = col_fun, show_row_dend = FALSE, show_column_dend = FALSE, 
		row_split = km, column_split = km,
		show_row_names = FALSE, show_column_names = FALSE,
		border = "#404040", row_title = NULL, 
		column_title = gt_render("B) Similarity matrix <span style='font-family:Times;'>**M**<sub>*b*</sub></span>"),
		row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_heatmap_legend = TRUE))
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>11</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 1)
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>12</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 2)
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>21</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 1)
	decorate_heatmap_body("Similarity", { ComplexHeatmap:::grid.text(gt_render("**M**<sub>22</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 2)
})


a = 0.95
pdf("figure1.pdf", width = 6.5, height = 3)
grid.newpage()
w = unit(1, "npc")
pushViewport(viewport(x = 0.25*w, width = 0.5*a*w))
grid.draw(plot1)
popViewport()
pushViewport(viewport(x = 0.75*w, width = 0.5*a*w))
grid.draw(plot2)
popViewport()
dev.off()

