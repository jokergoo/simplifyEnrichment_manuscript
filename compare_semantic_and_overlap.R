
library(GetoptLong)
library(ggplot2)
library(ComplexHeatmap)
library(simplifyEnrichment)
library(circlize)
library(cowplot)

setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity")

con_mat_list = list()
avg_n_list = list()
n_set_list = list()
for(i in 1:100) {
	qqcat("@{i}/100...\n")
	lt1 = readRDS(qq("../examples/random_BP/rds/clt_random_BP_@{i}.rds"))
	lt2 = readRDS(qq("../examples/random_BP_jaccard/rds/clt_random_BP_jaccard_@{i}.rds"))
	lt3 = readRDS(qq("../examples/random_BP_dice/rds/clt_random_BP_dice_@{i}.rds"))
	lt4 = readRDS(qq("../examples/random_BP_overlap/rds/clt_random_BP_overlap_@{i}.rds"))
	lt5 = readRDS(qq("../examples/random_BP_kappa/rds/clt_random_BP_kappa_@{i}.rds"))
	
	mat1 = lt1[[1]]
	mat2 = lt2[[1]]
	mat3 = lt3[[1]]
	mat4 = lt4[[1]]
	mat5 = lt5[[1]]

	df = rbind(data.frame(x = mat1[upper.tri(mat1)], cate = "semantic"),
		       data.frame(x = mat2[upper.tri(mat2)], cate = "jaccard"),
		       data.frame(x = mat3[upper.tri(mat3)], cate = "dice"),
		       data.frame(x = mat4[upper.tri(mat4)], cate = "overlap"),
		       data.frame(x = mat5[upper.tri(mat5)], cate = "kappa"))

	png(qq("random_BP/random_BP_@{i}_density_plot.png"), width = 1200*2, height = 600*2, res = 72*2)
	p1 = ggplot(df, aes(x=x, fill=factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")))) + 
		geom_density() + 
  		labs(title="Density Plot", x="Similarity") +
  		ggtitle(qq("random_BP_@{i}")) +
  		theme(legend.position = "none") +
  		facet_grid(~cate)
  	p2 = ggplot(df, aes(x=x, fill=factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")))) + 
		geom_density() + xlim(0, 0.05) +
  		labs(title="Density Plot", x="Similarity") +
  		ggtitle(qq("random_BP_@{i}, xlim = [0, 0.05]")) +
  		theme(legend.position = "none") +
  		facet_grid(~cate)
  	print(plot_grid(p1, p2, nrow = 2))
  	dev.off()

  	clt = list(semantic = as.character(lt1[[2]]$binary_cut),
  		       jaccard = as.character(lt2[[2]]$binary_cut),
  		       dice = as.character(lt3[[2]]$binary_cut),
  		       overlap = as.character(lt4[[2]]$binary_cut),
  		       kappa = as.character(lt5[[2]]$binary_cut)
  		   )

  	l = sapply(clt, function(x) any(is.na(x)))
	clt2 = clt[!l]
	
	ref_class = clt2$semantic

	clt2 = lapply(clt2, function(x) simplifyEnrichment:::relabel_class(x, ref_class, return_map = FALSE))
	clt2 = as.data.frame(clt2)

	png(qq("random_BP/random_BP_@{i}_concordance.png"), width = 150*2, height = 500*2, res = 72*2)
	ht = Heatmap(as.matrix(clt2), show_heatmap_legend = FALSE,
			row_order = do.call(order, clt), 
			width = unit(5, "mm")*ncol(clt2), column_names_rot = 45,
			row_title = qq("binary cut partitions from various similarity matrices, random_BP_@{i}"))
	draw(ht)
	dev.off()

	ht1 = grid.grabExpr(ht_clusters(mat1, clt[[1]], draw_word_cloud = FALSE, column_title = "semantic similarity"))
	ht2 = grid.grabExpr(ht_clusters(mat2, clt[[2]], draw_word_cloud = FALSE, column_title = "jaccard similarity"))
	ht3 = grid.grabExpr(ht_clusters(mat3, clt[[3]], draw_word_cloud = FALSE, column_title = "dice similarity"))
	ht4 = grid.grabExpr(ht_clusters(mat4, clt[[4]], draw_word_cloud = FALSE, column_title = "overlap similarity"))
	ht5 = grid.grabExpr(ht_clusters(mat5, clt[[5]], draw_word_cloud = FALSE, column_title = "kappa similarity"))

	jpeg(qq("random_BP/random_BP_@{i}_heatmap.jpg"), width = 1000, height = 600)
	print(plot_grid(ht1, ht2, ht3, ht4, ht5, nrow = 2))
	dev.off()

	html_file = qq("random_BP/random_BP_@{i}.html")
		link = qq("<a href='random_BP_@{i-1}.html' title='random_BP_@{i-1}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='random_BP_@{i+1}.html' title='random_BP_@{i+1}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/compare_similarity/random_BP_compare_similarity.html'>Back</a>")
		if(i == 1) {
			link = qq("<a href='random_BP_100.html' title='random_BP_100'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='random_BP_@{i+1}.html' title='random_BP_@{i+1}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/compare_similarity/random_BP_compare_similarity.html'>Back</a>")
		} else if(i == 100) {
			link = qq("<a href='random_BP_@{i-1}.html' title='random_BP_@{i-1}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='random_BP_1.html' title='random_BP_1'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/compare_similarity/random_BP_compare_similarity.html'>Back</a>")
		}
		writeLines(qq("
	<html>
	<head>
	<title>random_BP_@{i}</title>
	<link href='../main.css' rel='stylesheet'>
	</head>
	<body>
	<h2>random_BP_@{i}</h2>
	<hr />
	<p>@{link}</p>
	<p><b>Figure 1.</b> Binary clustering on similarity matrices from different measurements.</p>
	<p><img src='random_BP_@{i}_heatmap.jpg' /></p>
	<br>
	<p><b>Figure 2.</b> Compare clusterings.</p>
	<p><img src='random_BP_@{i}_concordance.png' width='200' /></p>
	<br>
	<p><b>Figure 3.</b> Density distribution of similarty values.</p>
	<p><img src='random_BP_@{i}_density_plot.png' width='1000' /></p>
	</body>
	</html>
	"), con = html_file)

	con_mat = simplifyEnrichment:::cmp_calc_concordance(clt2)
	avg_n = sapply(clt2, function(x) length(x)/length(unique(x)))
	n_set = sapply(clt2, function(x) length(unique(x)))

	con_mat_list[[i]] = con_mat
	avg_n_list[[i]] = avg_n
	n_set_list[[i]] = n_set
}

save(con_mat_list, avg_n_list, n_set_list, file = "rds/compare_similarity_random_BP.RData")

df1 = do.call(rbind, lapply(avg_n_list, function(x) data.frame(x = x, cate = names(x))))
df2 = do.call(rbind, lapply(n_set_list, function(x) data.frame(x = x, cate = names(x))))

png(qq("random_BP_compare_similarity.png"), width = 900*2, height = 300*2, res = 72*2)
p1 = ggplot(df1, aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x)) + 
	geom_boxplot() + labs(x = "", y = "Average number of terms per cluster") +
	scale_y_log10()
p2 = ggplot(df2, aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x)) + 
	geom_boxplot() + labs(x = "", y = "Number of clusters")


arr = array(dim = c(dim(con_mat), length(con_mat_list)))
for(i in seq_along(con_mat_list)) {
	arr[, , i] = con_mat_list[[i]]
}

avg_con_mat = apply(arr, 1:2, mean)
dimnames(avg_con_mat) = dimnames(con_mat)

ht = Heatmap(avg_con_mat, name = "Concordance", col = colorRamp2(c(0, 1), c("white", "red")),
	cell_fun = function(j, i, x, y, w, h, f) {
		grid.text(sprintf("%.3f", avg_con_mat[i, j]), x, y, gp = gpar(fontsize = 10))
	})
p3 = grid.grabExpr(draw(ht))

plot_grid(p2, p1, p3, nrow = 1 ,rel_widths = c(1, 1, 1.3))
dev.off()


tb = matrix(qq("<a href='/cmp_sim_random_BP/random_BP_@{1:100}.html'>random_BP_@{1:100}</a>", collapse = FALSE), 
	ncol = 10, byrow = TRUE)
writeLines(qq("
<html>
<head>
<title>Compare similarity measurements - random_BP</title>
<link href='main.css' rel='stylesheet'>
</head>
<body>
<h1>Compare similarity measurements - random_BP</h1>
<hr />
<p><b>Figure 1.</b> Left) Numbers of gene sets. Middle) Average number of terms per cluster. Right) Concordance of clusterings from different similarity matrices. The definition of the concordance score can be found <a href='/concordance.html'>here</a>.</p>
<p><img src='random_BP_compare_similarity.png' width='900' /></p>
<br>
<p><b>Table 1.</b> Details on individual datasets.</p>
@{knitr::kable(tb, row.names = FALSE, format = 'html', escape = FALSE)}
</body>
</html>
"), con = "random_BP_compare_similarity.html")


# =================================================

lt2 = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/rds/lt2_sim_all.rds")

compare = function(ont) {
	sim = lt2[[qq("@{ont}_sim_jaccard")]]
	n = length(sim)
	nm = names(sim)

	con_mat_list = list()
	avg_n_list = list()
	n_set_list = list()
	for(i in seq_len(n)) {
		qqcat("@{i}/@{n}\n")
		if(grepl("GO|DO", ont)) lt1 = readRDS(qq("../examples/EBI_Expression_Atlas_@{ont}/rds/clt_@{nm[i]}_@{ont}.rds"))
		lt2 = readRDS(qq("../examples/EBI_Expression_Atlas_@{ont}_jaccard/rds/clt_@{nm[i]}_@{ont}_jaccard.rds"))
		lt3 = readRDS(qq("../examples/EBI_Expression_Atlas_@{ont}_dice/rds/clt_@{nm[i]}_@{ont}_dice.rds"))
		lt4 = readRDS(qq("../examples/EBI_Expression_Atlas_@{ont}_overlap/rds/clt_@{nm[i]}_@{ont}_overlap.rds"))
		lt5 = readRDS(qq("../examples/EBI_Expression_Atlas_@{ont}_kappa/rds/clt_@{nm[i]}_@{ont}_kappa.rds"))
		
		if(grepl("GO|DO", ont)) mat1 = lt1[[1]]
		mat2 = lt2[[1]]
		mat3 = lt3[[1]]
		mat4 = lt4[[1]]
		mat5 = lt5[[1]]

		if(grepl("GO|DO", ont)) {
			df = rbind(data.frame(x = mat1[upper.tri(mat1)], cate = "semantic"),
			       data.frame(x = mat2[upper.tri(mat2)], cate = "jaccard"),
			       data.frame(x = mat3[upper.tri(mat3)], cate = "dice"),
			       data.frame(x = mat4[upper.tri(mat4)], cate = "overlap"),
			       data.frame(x = mat5[upper.tri(mat5)], cate = "kappa"))
		} else {
			df = rbind(data.frame(x = mat2[upper.tri(mat2)], cate = "jaccard"),
			       data.frame(x = mat3[upper.tri(mat3)], cate = "dice"),
			       data.frame(x = mat4[upper.tri(mat4)], cate = "overlap"),
			       data.frame(x = mat5[upper.tri(mat5)], cate = "kappa"))
		}

		dir.create(qq("EBI_Expression_Atlas_@{ont}"), showWarnings = FALSE)

		png(qq("EBI_Expression_Atlas_@{ont}/EBI_Expression_Atlas_@{ont}_@{nm[i]}_density_plot.png"), width = 1200*2, height = 600*2, res = 72*2)
		p1 = ggplot(df, aes(x=x, fill=factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")))) + 
			geom_density() + 
	  		labs(title="Density Plot", x="Similarity") +
	  		ggtitle(qq("EBI_Expression_Atlas_@{ont}_@{nm[i]}")) +
	  		theme(legend.position = "none") +
	  		facet_grid(~cate)
	  	p2 = ggplot(df, aes(x=x, fill=factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")))) + 
			geom_density() + xlim(0, 0.25) +
	  		labs(title="Density Plot", x="Similarity") +
	  		ggtitle(qq("EBI_Expression_Atlas_@{ont}_@{nm[i]}, xlim = [0, 0.25]")) +
	  		theme(legend.position = "none") +
	  		facet_grid(~cate)
	  	print(plot_grid(p1, p2, nrow = 2))
	  	dev.off()

	  	if(grepl("GO|DO", ont)) {
		  	clt = list(semantic = as.character(lt1[[2]]$binary_cut),
		  		       jaccard = as.character(lt2[[2]]$binary_cut),
		  		       dice = as.character(lt3[[2]]$binary_cut),
		  		       overlap = as.character(lt4[[2]]$binary_cut),
		  		       kappa = as.character(lt5[[2]]$binary_cut)
		  		   )
		} else {
			clt = list(jaccard = as.character(lt2[[2]]$binary_cut),
	  		       dice = as.character(lt3[[2]]$binary_cut),
	  		       overlap = as.character(lt4[[2]]$binary_cut),
	  		       kappa = as.character(lt5[[2]]$binary_cut)
	  		   )
		}

	  	l = sapply(clt, function(x) any(is.na(x)))
		clt2 = clt[!l]
		
		if(grepl("GO|DO", ont)) {
			ref_class = clt2$semantic
		} else {
			ref_class = clt2$jaccard
		}
		clt2 = lapply(clt2, function(x) simplifyEnrichment:::relabel_class(x, ref_class, return_map = FALSE))
		clt2 = as.data.frame(clt2)

		png(qq("EBI_Expression_Atlas_@{ont}/EBI_Expression_Atlas_@{ont}_@{nm[i]}_concordance.png"), width = 150*2, height = 500*2, res = 72*2)
		ht = Heatmap(as.matrix(clt2), show_heatmap_legend = FALSE,
				row_order = do.call(order, clt), 
				width = unit(5, "mm")*ncol(clt2), column_names_rot = 45)
		draw(ht)
		dev.off()

		if(grepl("GO|DO", ont)) {
			ht1 = grid.grabExpr(ht_clusters(mat1, clt[[1]], draw_word_cloud = FALSE, column_title = "semantic similarity"))
			ht2 = grid.grabExpr(ht_clusters(mat2, clt[[2]], draw_word_cloud = FALSE, column_title = "jaccard similarity"))
			ht3 = grid.grabExpr(ht_clusters(mat3, clt[[3]], draw_word_cloud = FALSE, column_title = "dice similarity"))
			ht4 = grid.grabExpr(ht_clusters(mat4, clt[[4]], draw_word_cloud = FALSE, column_title = "overlap similarity"))
			ht5 = grid.grabExpr(ht_clusters(mat5, clt[[5]], draw_word_cloud = FALSE, column_title = "kappa similarity"))
		} else {
			ht2 = grid.grabExpr(ht_clusters(mat2, clt[[1]], draw_word_cloud = FALSE, column_title = "jaccard similarity"))
			ht3 = grid.grabExpr(ht_clusters(mat3, clt[[2]], draw_word_cloud = FALSE, column_title = "dice similarity"))
			ht4 = grid.grabExpr(ht_clusters(mat4, clt[[3]], draw_word_cloud = FALSE, column_title = "overlap similarity"))
			ht5 = grid.grabExpr(ht_clusters(mat5, clt[[4]], draw_word_cloud = FALSE, column_title = "kappa similarity"))
		}

		
		if(grepl("GO|DO", ont)) {
			jpeg(qq("EBI_Expression_Atlas_@{ont}/EBI_Expression_Atlas_@{ont}_@{nm[i]}_heatmap.jpg"), width = 1000, height = 600)
			print(plot_grid(ht1, ht2, ht3, ht4, ht5, nrow = 2))
		} else {
			jpeg(qq("EBI_Expression_Atlas_@{ont}/EBI_Expression_Atlas_@{ont}_@{nm[i]}_heatmap.jpg"), width = 1200, height = 300)
			print(plot_grid(ht2, ht3, ht4, ht5, nrow = 1))
		}
		dev.off()

		html_file = qq("EBI_Expression_Atlas_@{ont}/EBI_Expression_Atlas_@{ont}_@{nm[i]}.html")
			link = qq("<a href='EBI_Expression_Atlas_@{ont}_@{nm[i-1]}.html' title='EBI_Expression_Atlas_@{ont}_@{nm[i-1]}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='EBI_Expression_Atlas_@{ont}_@{nm[i+1]}.html' title='EBI_Expression_Atlas_@{ont}_@{nm[i+1]}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/compare_similarity/EBI_Expression_Atlas_@{ont}_compare_similarity.html'>Back</a>")
			if(i == 1) {
				link = qq("<a href='EBI_Expression_Atlas_@{ont}_@{nm[n]}.html' title='EBI_Expression_Atlas_@{ont}_@{nm[n]}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='EBI_Expression_Atlas_@{ont}_@{nm[i+1]}.html' title='EBI_Expression_Atlas_@{ont}_@{nm[i+1]}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/compare_similarity/EBI_Expression_Atlas_@{ont}_compare_similarity.html'>Back</a>")
			} else if(i == n) {
				link = qq("<a href='EBI_Expression_Atlas_@{ont}_@{nm[i-1]}.html' title='EBI_Expression_Atlas_@{ont}_@{nm[i-1]}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='EBI_Expression_Atlas_@{ont}_@{nm[1]}.html' title='EBI_Expression_Atlas_@{ont}_@{nm[1]}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/compare_similarity/EBI_Expression_Atlas_@{ont}_compare_similarity.html'>Back</a>")
			}
			writeLines(qq("
		<html>
		<head>
		<title>EBI_Expression_Atlas_@{ont}_@{nm[i]}}</title>
		<link href='../main.css' rel='stylesheet'>
		</head>
		<body>
		<h2>EBI_Expression_Atlas_@{ont}_@{nm[i]}</h2>
		<hr />
		<p>@{link}</p>
		<p><b>Figure 1.</b> Binary clustering on similarity matrices from different measurements.</p>
		<p><img src='EBI_Expression_Atlas_@{ont}_@{nm[i]}_heatmap.jpg' /></p>
		<br>
		<p><b>Figure 2.</b> Compare clusterings.</p>
		<p><img src='EBI_Expression_Atlas_@{ont}_@{nm[i]}_concordance.png'  width='200' /></p>
		<br>
		<p><b>Figure 3.</b> Density distribution of similarty values.</p>
		<p><img src='EBI_Expression_Atlas_@{ont}_@{nm[i]}_density_plot.png' width='1000' /></p>
		</body>
		</html>
		"), con = html_file)

		con_mat = simplifyEnrichment:::cmp_calc_concordance(clt2)
		avg_n = sapply(clt2, function(x) length(x)/length(unique(x)))
		n_set = sapply(clt2, function(x) length(unique(x)))

		con_mat_list[[i]] = con_mat
		avg_n_list[[i]] = avg_n
		n_set_list[[i]] = n_set
	}

	save(con_mat_list, avg_n_list, n_set_list, file = qq("rds/compare_similarity_EBI_Expression_Atlas_@{ont}.RData"))

	df1 = do.call(rbind, lapply(avg_n_list, function(x) data.frame(x = x, cate = names(x))))
	df2 = do.call(rbind, lapply(n_set_list, function(x) data.frame(x = x, cate = names(x))))

	png(qq("EBI_Expression_Atlas_@{ont}_compare_similarity.png"), width = 900*2, height = 300*2, res = 72*2)
	p1 = ggplot(df1, aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x)) + 
		geom_boxplot() + labs(x = "", y = "Average number of terms per cluster") +
		scale_y_log10()
	p2 = ggplot(df2, aes(x = factor(cate, levels = c("semantic", "jaccard", "dice", "overlap", "kappa")), y = x)) + 
		geom_boxplot() + labs(x = "", y = "Number of clusters")

	arr = array(dim = c(dim(con_mat), length(con_mat_list)))
	for(i in seq_along(con_mat_list)) {
		arr[, , i] = con_mat_list[[i]]
	}

	avg_con_mat = apply(arr, 1:2, mean)
	dimnames(avg_con_mat) = dimnames(con_mat)

	ht = Heatmap(avg_con_mat, name = "concordance", col = colorRamp2(c(0, 1), c("white", "red")),
		cell_fun = function(j, i, x, y, w, h, f) {
			grid.text(sprintf("%.3f", avg_con_mat[i, j]), x, y, gp = gpar(fontsize = 10))
		})
	p3 = grid.grabExpr(draw(ht))

	print(plot_grid(p2, p1, p3, nrow = 1 ,rel_widths = c(1, 1, 1.3)))
	dev.off()



	tb = matrix(c(qq("<a href='/cmp_sim_EBI_Expression_Atlas_@{ont}/EBI_Expression_Atlas_@{ont}_@{nm}.html'>@{nm}</a>", collapse = FALSE), rep("", ceiling(length(nm)/5)*5 - length(nm))), 
		ncol = 5, byrow = TRUE)
	writeLines(qq("
	<html>
	<head>
	<title>Compare similarity measurements - EBI_Expression_Atlas_@{ont}</title>
	<link href='main.css' rel='stylesheet'>
	</head>
	<body>
	<h1>Compare similarity measurements - EBI_Expression_Atlas_@{ont}</h1>
	<hr />
	<p><b>Figure 1.</b> Left) Numbers of gene sets. Middle) Average number of terms per cluster. Right) Concordance of clusterings from different similarity matrices. The definition of the concordance score can be found <a href='/concordance.html'>here</a>.</p>
	<p><img src='EBI_Expression_Atlas_@{ont}_compare_similarity.png' width='900' /></p>
	<br>
	<p><b>Table 1.</b> Details on individual datasets.</p>
	@{knitr::kable(tb, row.names = FALSE, format = 'html', escape = FALSE)}
	</body>
	</html>
	"), con = qq("EBI_Expression_Atlas_@{ont}_compare_similarity.html"))

}

library(bsub)
bsub_opt$temp_dir = "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/bsub_temp"

for(ont in c("GO_BP", "DO", "KEGG", "Reactome", "MsigDB_C2_CGP", "MsigDB_C3_GTRD", "MsigDB_C3_MIR_Legacy", 
	"MsigDB_C3_MIRDB", "MsigDB_C3_TFT_Legacy", "MsigDB_C4_CGN", "MsigDB_C4_CM", "MsigDB_C7_")) {
	bsub_chunk(name = qq("cmp_sim_@{ont}"), hour = 10, memory = 20, variables = c("ont", "compare"),
	{	
		library(GetoptLong)
		library(ggplot2)
		library(ComplexHeatmap)
		library(simplifyEnrichment)
		library(circlize)
		library(cowplot)

		setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/compare_similarity")
		lt2 = readRDS("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/rds/lt2_sim_all.rds")
		compare(ont)
	})
}

# compare("GO_BP")
# compare("DO")
# compare("KEGG")
# compare("Reactome")
# compare("MsigDB_C2_CGP")
# compare("MsigDB_C3_GTRD")
# compare("MsigDB_C3_MIR_Legacy")
# compare("MsigDB_C3_MIRDB")
# compare("MsigDB_C3_TFT_Legacy")
# compare("MsigDB_C4_CGN")
# compare("MsigDB_C4_CM")
# compare("MsigDB_C7")

#######
writeLines(qq("
<html>
<head>
<title>Compare similarity measurements</title>
<link href='main.css' rel='stylesheet'>
</head>
<body>
<h1>Compare similarity measurements</h1>
<hr />
<ul>
<li><a href='random_BP_compare_similarity.html'>Compare similarity measurements - random_GO_BP</a></li>
<li><a href='EBI_Expression_Atlas_GO_BP_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_GO_BP</a></li>
<li><a href='EBI_Expression_Atlas_DO_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_DO</a></li>
<li><a href='EBI_Expression_Atlas_KEGG_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_KEGG</a></li>
<li><a href='EBI_Expression_Atlas_Reactome_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_Reactome</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C2_CGP_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C2_CGP</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C3_GTRD_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C3_GTRD</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C3_MIR_Legacy_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C3_MIR_Legacy</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C3_MIRDB_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C3_MIRDB</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C3_TFT_Legacy_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C3_TFT_Legacy</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C4_CGN_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C4_CGN</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C4_CM_compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C4_CM</a></li>
<li><a href='EBI_Expression_Atlas_MsigDB_C7__compare_similarity.html'>Compare similarity measurements - EBI_Expression_Atlas_MsigDB_C7</a></li>
</ul>
</body>
</html>
"), con = "compare_similarity.html")


servr::httd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/")

