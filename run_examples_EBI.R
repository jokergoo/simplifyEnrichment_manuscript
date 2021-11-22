
library(GetoptLong)

setwd("/omics/groups/OE0246/internal/guz/simplifyGO_test/examples")

# library(devtools)
# install_github("jokergoo/simplifyEnrichment")

library(simplifyEnrichment)
library(bsub)
bsub_opt$temp_dir = "/omics/groups/OE0246/internal/guz/simplifyGO_test/bsub_temp"

###################  EBI_Expression_Atlas  ########################

submit_EBI_Expression_Atlas = function(sim, onto) {
	n = length(sim)
	nm = names(sim)

	dir.create(qq("EBI_Expression_Atlas_@{onto}"), showWarnings = FALSE)
	dir.create(qq("EBI_Expression_Atlas_@{onto}/rds"), showWarnings = FALSE)
	dir.create(qq("EBI_Expression_Atlas_@{onto}/image"), showWarnings = FALSE)
	for(i in seq_len(n)) {
		# if(i > 10) break
		mat = as.matrix(sim[[i]])
		bsub_chunk(name = qq("Expression_Atlas_@{onto}_@{i}_@{nm[i]}"), variables = c("i", "mat", "n", "nm", "onto"), hour = ifelse(nrow(mat) < 500, 0.8, 4), memory = 10,
		{	
			library(GetoptLong)

   			setwd("/omics/groups/OE0246/internal/guz/simplifyGO_test/examples")

			library(simplifyEnrichment)
			
			qqcat("====== EBI_Expression_Atlas_@{onto} @{i}/@{n} ========\n")

			if(onto == "DO") {
				register_clustering_methods(
					binary_cut = function(mat, ...) binary_cut(mat, cutoff = 0.75, ...)
				)
			}
			
			clt = cmp_make_clusters(mat, method = all_clustering_methods())

			dir.create(qq("EBI_Expression_Atlas_@{onto}/image/"), showWarnings = FALSE, recursive = TRUE)
			dir.create(qq("EBI_Expression_Atlas_@{onto}/rds/"), showWarnings = FALSE, recursive = TRUE)

			saveRDS(list(mat, clt), file = qq("EBI_Expression_Atlas_@{onto}/rds/clt_@{nm[i]}_@{onto}.rds"))
			qqcat("rds saved to @{getwd()}/EBI_Expression_Atlas_@{onto}/rds/clt_@{nm[i]}_@{onto}.rds\n")

			jpeg(qq("EBI_Expression_Atlas_@{onto}/image/term_heatmap_@{nm[i]}_@{onto}.jpg"), width = 1800, height = 1800/5*2)
			cmp_make_plot(mat, clt, plot_type = "heatmap", nrow = 2)
			dev.off()

			png(qq("EBI_Expression_Atlas_@{onto}/image/term_heatmap_@{nm[i]}_@{onto}_barplot.png"), width = 800*2, height = 600*2, res = 72*2)
			cmp_make_plot(mat, clt)
			dev.off()

			png(qq("EBI_Expression_Atlas_@{onto}/image/term_heatmap_@{nm[i]}_@{onto}_binary_cut_select_cutoff.png"), width = 650*2, height = 600*2, res = 72*2)
			select_cutoff(mat)
			dev.off()

			cat("done.\n")
		}, output_dir = bsub_opt$output_dir, temp_dir = bsub_opt$temp_dir)
	}
}

bsub_opt$enforce = FALSE

lt2 = readRDS("/omics/groups/OE0246/internal/guz/simplifyGO_test/rds/lt2_sim_all.rds")

sapply(lt2, length)

submit_EBI_Expression_Atlas(lt2$GO_BP_sim, "GO_BP")
submit_EBI_Expression_Atlas(lt2$GO_BP_sim_kappa, "GO_BP_kappa")
submit_EBI_Expression_Atlas(lt2$GO_BP_sim_jaccard, "GO_BP_jaccard")
submit_EBI_Expression_Atlas(lt2$GO_BP_sim_dice, "GO_BP_dice")
submit_EBI_Expression_Atlas(lt2$GO_BP_sim_overlap, "GO_BP_overlap")
submit_EBI_Expression_Atlas(lt2$DO_sim, "DO")
submit_EBI_Expression_Atlas(lt2$DO_sim_kappa, "DO_kappa")
submit_EBI_Expression_Atlas(lt2$DO_sim_jaccard, "DO_jaccard")
submit_EBI_Expression_Atlas(lt2$DO_sim_dice, "DO_dice")
submit_EBI_Expression_Atlas(lt2$DO_sim_overlap, "DO_overlap")
submit_EBI_Expression_Atlas(lt2$KEGG_sim_kappa, "KEGG_kappa")
submit_EBI_Expression_Atlas(lt2$KEGG_sim_jaccard, "KEGG_jaccard")
submit_EBI_Expression_Atlas(lt2$KEGG_sim_dice, "KEGG_dice")
submit_EBI_Expression_Atlas(lt2$KEGG_sim_overlap, "KEGG_overlap")
submit_EBI_Expression_Atlas(lt2$Reactome_sim_kappa, "Reactome_kappa")
submit_EBI_Expression_Atlas(lt2$Reactome_sim_jaccard, "Reactome_jaccard")
submit_EBI_Expression_Atlas(lt2$Reactome_sim_dice, "Reactome_dice")
submit_EBI_Expression_Atlas(lt2$Reactome_sim_overlap, "Reactome_overlap")

for(sm in c("C2_CGP", "C3_MIR_Legacy", "C3_TFT_Legacy", "C3_GTRD", "C3_MIRDB", "C4_CM", "C4_CGN", "C7_IMMUNESIGDB")) {
	submit_EBI_Expression_Atlas(lt2[[qq("MsigDB_@{sm}_sim_kappa")]], qq("MsigDB_@{sm}_kappa"))
	submit_EBI_Expression_Atlas(lt2[[qq("MsigDB_@{sm}_sim_jaccard")]], qq("MsigDB_@{sm}_jaccard"))
	submit_EBI_Expression_Atlas(lt2[[qq("MsigDB_@{sm}_sim_dice")]], qq("MsigDB_@{sm}_dice"))
	submit_EBI_Expression_Atlas(lt2[[qq("MsigDB_@{sm}_sim_overlap")]], qq("MsigDB_@{sm}_overlap"))
}

generate_EBI_Expression_Atlas_html = function(sim, onto) {
	n = length(sim)
	nm = names(sim)

	info = ifelse(!grepl('jaccard|dice|overlap|kappa', onto), "applied on semantic similarity matrix", 
		qq("applied on @{gsub('^.*(jaccard|dice|overlap|kappa).*$', '\\\\1', onto)} coeffcient matrix"))

	file.copy("main.css", qq("EBI_Expression_Atlas_@{onto}/main.css"))

	m = NULL
	for(i in seq_len(n)) {
		qqcat("reading EBI_Expression_Atlas_@{onto}/rds/clt_@{nm[i]}_@{onto}.rds, @{i}/@{n}\n")
		oe = try(lt <- readRDS(qq("EBI_Expression_Atlas_@{onto}/rds/clt_@{nm[i]}_@{onto}.rds")))

		if(inherits(oe, "try-error")) {
			m = rbind(m, rep(NA, 11))
		} else {
			m = rbind(m, sapply(lt[[2]], function(x) {
				if(identical(x, NA)) {
					"NA"
				} else {
					tb = table(x)
					paste0(length(tb), "(", sum(tb >= 5), ")")
				}
			}))

			html_file = qq("EBI_Expression_Atlas_@{onto}/image/term_heatmap_@{nm[i]}_@{onto}.html")
			link = qq("<a href='term_heatmap_@{nm[i-1]}_@{onto}.html' title='term_heatmap_@{nm[i-1]}_@{onto}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='term_heatmap_@{nm[i+1]}_@{onto}.html' title='term_heatmap_@{nm[i+1]}_@{onto}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/examples/EBI_Expression_Atlas_@{onto}.html'>Back</a>")
			if(i == 1) {
				link = qq("<a href='term_heatmap_@{nm[n]}_@{onto}.html' title='term_heatmap_@{nm[n]}_@{onto}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='term_heatmap_@{nm[i+1]}_@{onto}.html' title='term_heatmap_@{nm[i+1]}_@{onto}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/examples/EBI_Expression_Atlas_@{onto}.html'>Back</a>")
			} else if(i == n) {
				link = qq("<a href='term_heatmap_@{nm[i-1]}_@{onto}.html' title='term_heatmap_@{nm[i-1]}_@{onto}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='term_heatmap_@{nm[1]}_@{onto}.html' title='term_heatmap_@{nm[1]}_@{onto}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/examples/EBI_Expression_Atlas_@{onto}.html'>Back</a>")
			}
			writeLines(qq("
<html>
<head>
<title>term_heatmap_@{nm[i]}_@{onto}</title>
<link href='../main.css' rel='stylesheet'>
</head>
<body>
<h2>term_heatmap_@{nm[i]}_@{onto} (@{info})</h2>
<hr />
<p>@{link}</p>
<p><b>Figure 1.</b> Clusterings from different methods. Green bars represent the small clusters with number of elements less than 5.</p>
<p><img src='term_heatmap_@{nm[i]}_@{onto}.jpg' width='1200' /></p>
<br>
<p><b>Figure 2.</b> Compare clustering results. Top left panel: Similarity heatmap with clusterings from different methods. Bottom left panel: Concordance between clustering methods. The concordance measures how similar two clusterings are. Right panel: The difference score, number of clusters and the block mean of different clusterings.</p>
<p><img src='term_heatmap_@{nm[i]}_@{onto}_barplot.png' width='800' /></p>
<br>
<p><b>Figure 3.</b> Select proper cutoff for binary cut.</p>
<p><img src='term_heatmap_@{nm[i]}_@{onto}_binary_cut_select_cutoff.png' width='600' /></p>
</body>
</html>
"), con = html_file)
		}
	}

	url = gsub('^.*(E-[A-Z]+-[0-9]+)_.*$', '\\1', nm)
	m = cbind(ID = qq("<a href='https://www.ebi.ac.uk/gxa/experiments/@{url}/Results'>@{nm}</a>", collapse = FALSE), m)
	m = cbind(m, "Details" = qq("<a href='https://simplifyEnrichment.github.io/EBI_Expression_Atlas_@{onto}/image/term_heatmap_@{nm}_@{onto}.html'>view</a>", collapse = FALSE))

	# the main page
	writeLines(qq("
	<html>
	<head>
	<title>EBI_Expression_Atlas_@{onto}</title>
	<link href='main.css' rel='stylesheet'>
	</head>
	<body>
	<h1>EBI_Expression_Atlas_@{onto} (@{info})</h1>
	<hr>
	<p><b>Figure 1.</b>Compare clustering results. Left panel: The difference score, number of clusters and the block mean of different clusterings. Right panel: Concordance between clustering methods. The concordance measures how similar two clusterings are. The definition of the concordance score can be found <a href='/concordance.html'>here</a>.</p>
	<p><img src='EBI_Expression_Atlas_@{onto}_results.png' width='800' /></p>
	<p><b>Table 1.</b>Number of clusters identified by each clustering method. Numbers in the table indicate the number of clusters. The numbers inside the parentheses are the number of clusters with size >= 5. </p>
	@{knitr::kable(m, row.names = FALSE, format = 'html', escape = FALSE)}
	</body>
	</html>
	"), con = qq("EBI_Expression_Atlas_@{onto}.html"))
}

generate_EBI_Expression_Atlas_html(lt2$GO_BP_sim, "GO_BP")
generate_EBI_Expression_Atlas_html(lt2$GO_BP_sim_kappa, "GO_BP_kappa")
generate_EBI_Expression_Atlas_html(lt2$GO_BP_sim_jaccard, "GO_BP_jaccard")
generate_EBI_Expression_Atlas_html(lt2$GO_BP_sim_dice, "GO_BP_dice")
generate_EBI_Expression_Atlas_html(lt2$GO_BP_sim_overlap, "GO_BP_overlap")
generate_EBI_Expression_Atlas_html(lt2$DO_sim, "DO")
generate_EBI_Expression_Atlas_html(lt2$DO_sim_kappa, "DO_kappa")
generate_EBI_Expression_Atlas_html(lt2$DO_sim_jaccard, "DO_jaccard")
generate_EBI_Expression_Atlas_html(lt2$DO_sim_dice, "DO_dice")
generate_EBI_Expression_Atlas_html(lt2$DO_sim_overlap, "DO_overlap")
generate_EBI_Expression_Atlas_html(lt2$KEGG_sim_kappa, "KEGG_kappa")
generate_EBI_Expression_Atlas_html(lt2$KEGG_sim_jaccard, "KEGG_jaccard")
generate_EBI_Expression_Atlas_html(lt2$KEGG_sim_dice, "KEGG_dice")
generate_EBI_Expression_Atlas_html(lt2$KEGG_sim_overlap, "KEGG_overlap")
generate_EBI_Expression_Atlas_html(lt2$Reactome_sim_kappa, "Reactome_kappa")
generate_EBI_Expression_Atlas_html(lt2$Reactome_sim_jaccard, "Reactome_jaccard")
generate_EBI_Expression_Atlas_html(lt2$Reactome_sim_dice, "Reactome_dice")
generate_EBI_Expression_Atlas_html(lt2$Reactome_sim_overlap, "Reactome_overlap")

for(sm in c("C2_CGP", "C3_MIR_Legacy", "C3_TFT_Legacy", "C3_GTRD", "C3_MIRDB", "C4_CM", "C4_CGN", "C7_IMMUNESIGDB")) {
	generate_EBI_Expression_Atlas_html(lt2[[qq("MsigDB_@{sm}_sim_kappa")]], qq("MsigDB_@{sm}_kappa"))
	generate_EBI_Expression_Atlas_html(lt2[[qq("MsigDB_@{sm}_sim_jaccard")]], qq("MsigDB_@{sm}_jaccard"))
	generate_EBI_Expression_Atlas_html(lt2[[qq("MsigDB_@{sm}_sim_dice")]], qq("MsigDB_@{sm}_dice"))
	generate_EBI_Expression_Atlas_html(lt2[[qq("MsigDB_@{sm}_sim_overlap")]], qq("MsigDB_@{sm}_overlap"))
}

servr::httd("/omics/groups/OE0246/internal/guz/simplifyGO_test/")

######################### summarize #####################

library(GetoptLong)
library(ggplot2)
library(grid)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

summarize = function(template, sim, output, rerun = FALSE) {

	if(!file.exists(output) || rerun) {
		if(is.atomic(sim)) {
			n = sim
		} else {
			n = length(sim)
			nm = names(sim)
		}

		concordance_mat_lt = list()
		diff_score_df = data.frame(method = character(0), value = numeric(0))
		cluster_number_df = data.frame(method = character(0), value = numeric(0))
		cluster_number2_df = data.frame(method = character(0), value = numeric(0))
		block_mean_df = data.frame(method = character(0), value = numeric(0))
		for(i in 1:n) {
			qqcat("reading @{qq(template)}, @{i}/@{n}\n")
			oe = try(lt <- readRDS(qq(template)))
			# if(inherits(oe, "try-error")) next

			concordance_mat_lt[[i]] = simplifyEnrichment:::cmp_calc_concordance(lt[[2]])
			
			x = sapply(lt[[2]], function(x) {
				if(is.na(x[1])) {
					NA
				} else if(length(unique(x)) == 1) {
					NA
				} else {
					difference_score(lt[[1]], x)
				}
			})
			diff_score_df = rbind(diff_score_df, data.frame(method = names(x), value = x))
			
			x = sapply(lt[[2]], function(x) {
				if(is.na(x[1])) {
					NA
				} else {
					length(table(x))
				}
			})
			cluster_number_df = rbind(cluster_number_df, data.frame(method = names(x), value = x))
			
			x = sapply(lt[[2]], function(x) {
				if(is.na(x[1])) {
					NA
				} else {
					tb = table(x); 
					sum(tb >= 5)
				}
			})
			cluster_number2_df = rbind(cluster_number2_df, data.frame(method = names(x), value = x))
			
			x = sapply(lt[[2]], function(x) {
				if(is.na(x[1])) {
					NA
				} else {
					simplifyEnrichment:::block_mean(lt[[1]], x)
				}
			})
			block_mean_df = rbind(block_mean_df, data.frame(method = names(x), value = x))
		}

		diff_score_df[, 1] = factor(diff_score_df[, 1], levels = all_clustering_methods())
		cluster_number_df[, 1] = factor(cluster_number_df[, 1], levels = all_clustering_methods())
		cluster_number2_df[, 1] = factor(cluster_number2_df[, 1], levels = all_clustering_methods())
		block_mean_df[, 1] = factor(block_mean_df[, 1], levels = all_clustering_methods())


		arr = array(dim = c(dim(concordance_mat_lt[[1]]), length(concordance_mat_lt)))
		dimnames(arr) = c(dimnames(concordance_mat_lt[[1]]), list(1:length(concordance_mat_lt)))
		for(i in 1:length(concordance_mat_lt)) {
			arr[, , i] = concordance_mat_lt[[i]]
		}

		saveRDS(list(diff_score_df = diff_score_df,
			         cluster_number_df = cluster_number_df,
			         cluster_number2_df = cluster_number2_df,
			         block_mean_df = block_mean_df,
			         arr = arr),
		        file = output,
		        compress = "xz")
	} else {
		lt_foo = readRDS(output)
		diff_score_df = lt_foo$diff_score_df
		cluster_number_df = lt_foo$cluster_number_df
		cluster_number2_df = lt_foo$cluster_number2_df
		block_mean_df = lt_foo$block_mean_df
		arr = lt_foo$arr
	}

	p1 = ggplot(diff_score_df, aes(x = method, y = value)) +
		geom_boxplot() + ylab("Difference score") +
		theme(axis.title.x = element_blank(), axis.text.x = element_blank())
	cluster_number_df$type = "All sizes"
	cluster_number2_df$type = "Size >= 5"
	p2 = ggplot(rbind(cluster_number_df, cluster_number2_df), aes(x = method, y = value, col = type)) +
		geom_boxplot() + ylab("Number of clusters") + labs(col = "Type") +
		theme(axis.title.x = element_blank(), axis.text.x = element_blank())
	p3 = ggplot(block_mean_df, aes(x = method, y = value)) +
		geom_boxplot() + ylab("Block mean") +
		theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

	m_mean = apply(arr, 1:2, mean, na.rm = TRUE)
	m_mean[is.na(m_mean)] = 0
	p4 = grid.grabExpr(draw(Heatmap(m_mean, name = "Concordance", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
		column_names_rot = 45)))

	png(gsub("rds$", "png", output), width = 800*2, height = 400*2, res = 72*2)
	print(plot_grid(plot_grid(p1, p2, p3, nrow = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 1.5)),
		p4, nrow = 1))
	dev.off()

}


summarize("EBI_Expression_Atlas_GO_BP/rds/clt_@{nm[i]}_GO_BP.rds", lt2$GO_BP_sim, output = "EBI_Expression_Atlas_GO_BP_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_GO_BP_kappa/rds/clt_@{nm[i]}_GO_BP_kappa.rds", lt2$GO_BP_sim_kappa, output = "EBI_Expression_Atlas_GO_BP_kappa_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_GO_BP_jaccard/rds/clt_@{nm[i]}_GO_BP_jaccard.rds", lt2$GO_BP_sim_jaccard, output = "EBI_Expression_Atlas_GO_BP_jaccard_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_GO_BP_dice/rds/clt_@{nm[i]}_GO_BP_dice.rds", lt2$GO_BP_sim_dice, output = "EBI_Expression_Atlas_GO_BP_dice_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_GO_BP_overlap/rds/clt_@{nm[i]}_GO_BP_overlap.rds", lt2$GO_BP_sim_overlap, output = "EBI_Expression_Atlas_GO_BP_overlap_results.rds", rerun = TRUE)

summarize("EBI_Expression_Atlas_DO/rds/clt_@{nm[i]}_DO.rds", lt2$DO_sim, output = "EBI_Expression_Atlas_DO_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_DO_kappa/rds/clt_@{nm[i]}_DO_kappa.rds", lt2$DO_sim_kappa, output = "EBI_Expression_Atlas_DO_kappa_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_DO_jaccard/rds/clt_@{nm[i]}_DO_jaccard.rds", lt2$DO_sim_jaccard, output = "EBI_Expression_Atlas_DO_jaccard_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_DO_dice/rds/clt_@{nm[i]}_DO_dice.rds", lt2$DO_sim_dice, output = "EBI_Expression_Atlas_DO_dice_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_DO_overlap/rds/clt_@{nm[i]}_DO_overlap.rds", lt2$DO_sim_overlap, output = "EBI_Expression_Atlas_DO_overlap_results.rds", rerun = TRUE)

summarize("EBI_Expression_Atlas_KEGG_kappa/rds/clt_@{nm[i]}_KEGG_kappa.rds", lt2$KEGG_sim_kappa, output = "EBI_Expression_Atlas_KEGG_kappa_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_KEGG_jaccard/rds/clt_@{nm[i]}_KEGG_jaccard.rds", lt2$KEGG_sim_jaccard, output = "EBI_Expression_Atlas_KEGG_jaccard_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_KEGG_dice/rds/clt_@{nm[i]}_KEGG_dice.rds", lt2$KEGG_sim_dice, output = "EBI_Expression_Atlas_KEGG_dice_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_KEGG_overlap/rds/clt_@{nm[i]}_KEGG_overlap.rds", lt2$KEGG_sim_overlap, output = "EBI_Expression_Atlas_KEGG_overlap_results.rds", rerun = TRUE)

summarize("EBI_Expression_Atlas_Reactome_kappa/rds/clt_@{nm[i]}_Reactome_kappa.rds", lt2$Reactome_sim_kappa, output = "EBI_Expression_Atlas_Reactome_kappa_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_Reactome_jaccard/rds/clt_@{nm[i]}_Reactome_jaccard.rds", lt2$Reactome_sim_jaccard, output = "EBI_Expression_Atlas_Reactome_jaccard_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_Reactome_dice/rds/clt_@{nm[i]}_Reactome_dice.rds", lt2$Reactome_sim_dice, output = "EBI_Expression_Atlas_Reactome_dice_results.rds", rerun = TRUE)
summarize("EBI_Expression_Atlas_Reactome_overlap/rds/clt_@{nm[i]}_Reactome_overlap.rds", lt2$Reactome_sim_overlap, output = "EBI_Expression_Atlas_Reactome_overlap_results.rds", rerun = TRUE)


for(sm in c("C2_CGP", "C3_MIR_Legacy", "C3_TFT_Legacy", "C3_GTRD", "C3_MIRDB", "C4_CM", "C4_CGN", "C7_IMMUNESIGDB")) {
	for(coef in c("kappa", "jaccard", "dice", "overlap")) {
		bsub_chunk(name = qq("summarize_@{sm}_@{coef}"), hour = 10, memory = 30, variables = c("sm", "coef", "summarize"),
		{	
			library(GetoptLong)
			library(ggplot2)
			library(grid)
			library(cowplot)
			library(ComplexHeatmap)
			library(circlize)
			library(simplifyEnrichment)

			setwd("/omics/groups/OE0246/internal/guz/simplifyGO_test/examples")
			lt2 = readRDS("/omics/groups/OE0246/internal/guz/simplifyGO_test/rds/lt2_sim_all.rds")
			
			summarize(qq("EBI_Expression_Atlas_MsigDB_`sm`_@{coef}/rds/clt_@{nm[i]}_MsigDB_`sm`_`coef`.rds", code.pattern = "`CODE`"), 
				lt2[[qq("MsigDB_@{sm}_sim_@{coef}")]], output = qq("EBI_Expression_Atlas_MsigDB_@{sm}_@{coef}_results.rds"), rerun = TRUE)
		})
	}
}

