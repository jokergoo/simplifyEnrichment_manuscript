setwd("~/manuscript/simplifyEnrichment")

library(simplifyEnrichment)

set.seed(888)
go_id = random_GO(500)
mat = GO_similarity(go_id)

clt = cmp_make_clusters(mat, method = all_clustering_methods())
clt = as.data.frame(clt)
methods = names(clt)


library(ComplexHeatmap)
library(GetoptLong)
nrow = 3

pdf("figure4.pdf", width = 18, height = 18/4*3)

pl = list()
lgd = NULL

for(i in seq_along(methods)) {
	if(any(is.na(clt[[i]]))) {
		pl[[i]] = textGrob(qq("@{methods[i]}\nan error occured."))
	} else {
		pl[[i]] = grid.grabExpr(ht <- ht_clusters(mat, clt[[i]], draw_word_cloud = FALSE, 
			column_title = qq("@{LETTERS[i]}) @{nrow(mat)} terms clustered by '@{methods[i]}'"),
			show_heatmap_legend = FALSE))
		lgd1 = color_mapping_legend(ht@ht_list[[1]]@matrix_color_mapping, plot = FALSE,
			legend_direction = "horizontal", title_position = "lefttop")
		lgd2 = Legend(labels = "All other small clusters (size < 5)", legend_gp = gpar(fill = "darkgreen"))
		lgd = packLegend(lgd1, lgd2, direction = "horizontal")
	}
}

n_all = sapply(clt, function(x) {
	if(any(is.na(x))) {
		NA
	} else {
		length(unique(x))
	}
})
n_big = sapply(clt, function(x) {
	if(any(is.na(x))) {
		NA
	} else {
		tb = table(x)
		sum(tb >= 5)
	}
})

tb = data.frame(Method = methods, "All clusters" = n_all, "Large clusters (size >= 5)" = n_big, check.names = FALSE)

pl[[length(pl) + 1]] = gridExtra::tableGrob(tb, rows = NULL)

np = length(pl)

ncol = ceiling(np/nrow)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol)))
for(i in 1:np) {
	ir = ceiling(i/ncol)
	ic = i %% ncol; if(ic == 0) ic = ncol
	pushViewport(viewport(layout.pos.row = ir, layout.pos.col = ic))
	if(i < np) {
		grid.draw(pl[[i]])
	} else {
		pushViewport(viewport(x = unit(0, "npc") + unit(2, "mm"), y = unit(1, "npc") - unit(1, "cm"),
			width = sum(pl[[i]]$widths), height = sum(pl[[i]]$heights),
			just = c("left", "top")))
		grid.draw(pl[[i]])
		grid.text("L) Numbers of clusters", y = unit(1, "npc") + unit(2.5, "mm"), just = "bottom", gp = gpar(fontsize = 14))
		popViewport()

		if(!is.null(lgd)) {
			draw(lgd, x = unit(0, "npc") + unit(2, "mm"), y = unit(1, "npc") - sum(pl[[i]]$heights) - unit(1.5, "cm"), just = c("left", "top"))
		}
	}
	popViewport()
}
	
dev.off()
