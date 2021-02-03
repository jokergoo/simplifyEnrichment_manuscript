library(simplifyEnrichment)

setwd("~/manuscript/simplifyEnrichment")

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
	draw(Heatmap(mat1, name = "mat1", col = col_fun, show_row_dend = FALSE, show_column_dend = FALSE, 
		row_split = km, column_split = km,
		show_row_names = FALSE, show_column_names = FALSE,
		border = "#404040", row_title = NULL, 
		column_title = gt_render("A) Similarity matrix <span style='font-family:Times;'>**M**<sub>*a*</sub></span>"),
		row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_heatmap_legend = FALSE))
	decorate_heatmap_body("mat1", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>11</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 1)
	decorate_heatmap_body("mat1", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>12</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 2)
	decorate_heatmap_body("mat1", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>21</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 1)
	decorate_heatmap_body("mat1", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>22</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 2)
})
mat2 = mat[ind2, ind2]
km = kmeans(mat2, centers = 2)$cluster
plot2 = grid.grabExpr({
	draw(Heatmap(mat2, name = "mat2", col = col_fun, show_row_dend = FALSE, show_column_dend = FALSE, 
		row_split = km, column_split = km,
		show_row_names = FALSE, show_column_names = FALSE,
		border = "#404040", row_title = NULL, 
		column_title = gt_render("B) Similarity matrix <span style='font-family:Times;'>**M**<sub>*b*</sub></span>"),
		row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), show_heatmap_legend = FALSE))
	decorate_heatmap_body("mat2", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>11</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 1)
	decorate_heatmap_body("mat2", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>12</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 1, column_slice = 2)
	decorate_heatmap_body("mat2", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>21</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 1)
	decorate_heatmap_body("mat2", { ComplexHeatmap:::grid.text(gt_render("*s*<sub>22</sub>"), x = 0.5, y = 0.5, gp = gpar(fontfamily = "Times", fontsize = 20)) }, row_slice = 2, column_slice = 2)
})


lgd = Legend(title = "Similarity", col_fun = col_fun)

a = 0.95
pdf("figure1.pdf", width = 6, height = 3)
grid.newpage()
w = unit(1, "npc") - grobWidth(lgd@grob) - unit(4, "mm")
pushViewport(viewport(x = 0.25*w, width = 0.5*a*w))
grid.draw(plot1)
popViewport()
pushViewport(viewport(x = 0.75*w, width = 0.5*a*w))
grid.draw(plot2)
popViewport()
draw(lgd, x = unit(1, "npc") - unit(2, "mm"), just = "right")
dev.off()

