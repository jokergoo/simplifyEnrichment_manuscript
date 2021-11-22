library(GetoptLong)

setwd("/omics/groups/OE0246/internal/guz/simplifyGO_test/examples")


library(simplifyEnrichment)
library(bsub)
bsub_opt$temp_dir = "/omics/groups/OE0246/internal/guz/simplifyGO_test/bsub_temp"

######################## random GO ###################

bsub_opt$enforce = TRUE
for(ont in c("BP", "MF", "CC")) {
	for(i in 1:100) {
	seed = round(runif(1)*2^30)
	bsub_chunk(name = qq("random_@{ont}_@{i}"), variables = c("i", "ont", "seed"), memory = 10, hour = ifelse(ont == "BP", 3, 0.5),
	{

		library(GetoptLong)

	    setwd("/omics/groups/OE0246/internal/guz/simplifyGO_test/examples")

		library(simplifyEnrichment)

		set.seed(seed)
		qqcat("random seed: @{seed}\n")

		qqcat("====== random GO_@{ont} @{i} ========\n")
		go_id = random_GO(500, ont = ont)
		mat = GO_similarity(go_id)

		clt = cmp_make_clusters(mat, method = all_clustering_methods())

		dir.create(qq("random_@{ont}/image/"), showWarnings = TRUE, recursive = TRUE)
		dir.create(qq("random_@{ont}/rds/"), showWarnings = TRUE, recursive = TRUE)

		jpeg(qq("random_@{ont}/image/GO_heatmap_random_@{ont}_@{i}.jpg"), width = 1800, height = 1800/5*2)
		cmp_make_plot(mat, clt, nrow = 2, plot_type = "heatmap")
		dev.off()

		png(qq("random_@{ont}/image/GO_heatmap_random_@{ont}_@{i}_barplot.png"), width = 800*2, height = 600*2, res = 72*2)
		cmp_make_plot(mat, clt)
		dev.off()

		png(qq("random_@{ont}/image/GO_heatmap_random_@{ont}_@{i}_binary_cut_select_cutoff.png"), width = 650*2, height = 600*2, res = 72*2)
		select_cutoff(mat)
		dev.off()

		saveRDS(list(mat, clt, seed), file = qq("random_@{ont}/rds/clt_random_@{ont}_@{i}.rds"))

		if(ont == "BP") {
			GO_DATA = clusterProfiler:::get_GO_data("org.Hs.eg.db", ont, "ENTREZID")
			gl = GO_DATA$PATHID2EXTID[which(names(GO_DATA$PATHID2EXTID) %in% go_id)]
			for(method in c("jaccard", "kappa", "dice", "overlap")) {
				mat = term_similarity(gl, method = method)
				clt = cmp_make_clusters(mat, method = all_clustering_methods())

				dir.create(qq("random_@{ont}_@{method}/image/"), showWarnings = TRUE, recursive = TRUE)
				dir.create(qq("random_@{ont}_@{method}/rds/"), showWarnings = TRUE, recursive = TRUE)

				jpeg(qq("random_@{ont}_@{method}/image/GO_heatmap_random_@{ont}_@{method}_@{i}.jpg"), width = 1800, height = 1800/5*2)
				cmp_make_plot(mat, clt, nrow = 2, plot_type = "heatmap")
				dev.off()

				png(qq("random_@{ont}_@{method}/image/GO_heatmap_random_@{ont}_@{method}_@{i}_barplot.png"), width = 800*2, height = 600*2, res = 72*2)
				cmp_make_plot(mat, clt)
				dev.off()

				png(qq("random_@{ont}_@{method}/image/GO_heatmap_random_@{ont}_@{method}_@{i}_binary_cut_select_cutoff.png"), width = 650*2, height = 600*2, res = 72*2)
				select_cutoff(mat)
				dev.off()

				saveRDS(list(mat, clt, seed), file = qq("random_@{ont}_@{method}/rds/clt_random_@{ont}_@{method}_@{i}.rds"))
			}
		}

		cat("done.\n")
	})
	}
}


for(ont in c("BP", "MF", "CC")) {
	for(i in 1:100) {
		if(!is_job_finished(qq("random_@{ont}_@{i}"))) {
			qqcat("not finished: random_@{ont}_@{i}\n")
		}
	}
}

for(ont in c("BP", "MF", "CC")) {
	for(i in 1:100) {
		status = job_status_by_name(qq("random_@{ont}_@{i}"))
		if(status != "DONE") qqcat("random_@{ont}_@{i}... @{status}\n")
	}
}

setwd("/omics/groups/OE0246/internal/guz/simplifyGO_test/examples")

## generate HTML page for every dataset
for(ont in c("BP", "MF", "CC")) {
	for(method in c("", "jaccard", "kappa", "dice", "overlap")) {

		if(ont %in% c("CC", "MF") && method != "") next

		if(method == "") {
			substr = ""
		} else {
			substr = paste0("_", method)
		}
		info = ifelse(substr == '', "applied on semantic similarity matrix", qq("applied on @{gsub('_', '', substr)} coeffcient matrix"))

		file.copy("main.css", qq("random_@{ont}@{substr}/main.css"))

		m = NULL
		for(i in 1:100) {
			qqcat("reading random_@{ont}@{substr}/rds/clt_random_@{ont}@{substr}_@{i}.rds, @{i}/100\n")
			lt = readRDS(qq("random_@{ont}@{substr}/rds/clt_random_@{ont}@{substr}_@{i}.rds"))

			html_file = qq("random_@{ont}@{substr}/image/GO_heatmap_random_@{ont}@{substr}_@{i}.html")
			link = qq("<a href='GO_heatmap_random_@{ont}@{substr}_@{i-1}.html' title='GO_heatmap_random_@{ont}@{substr}_@{i-1}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='GO_heatmap_random_@{ont}@{substr}_@{i+1}.html' title='GO_heatmap_random_@{ont}@{substr}_@{i+1}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/examples/random_@{ont}@{substr}.html'>Back</a>")
			if(i == 1) {
				link = qq("<a href='GO_heatmap_random_@{ont}@{substr}_100.html' title='GO_heatmap_random_@{ont}@{substr}_100'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='GO_heatmap_random_@{ont}@{substr}_@{i+1}.html' title='GO_heatmap_random_@{ont}@{substr}_@{i+1}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/examples/random_@{ont}@{substr}.html'>Back</a>")
			} else if(i == 100) {
				link = qq("<a href='GO_heatmap_random_@{ont}@{substr}_@{i-1}.html' title='GO_heatmap_random_@{ont}@{substr}_@{i-1}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='GO_heatmap_random_@{ont}@{substr}_1.html' title='GO_heatmap_random_@{ont}@{substr}_1'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='/examples/random_@{ont}@{substr}.html'>Back</a>")
			}
			writeLines(qq("
		<html>
		<head>
		<title>GO_heatmap_random_@{ont}@{substr}_@{i}</title>
		<link href='../main.css' rel='stylesheet'>
		</head>
		<body>
		<h2>GO_heatmap_random_@{ont}@{substr}_@{i} (@{info})</h2>
		<hr />
		<p>@{link}</p>
		<p><b>Figure 1.</b> Clusterings from different methods. Green bars represent the small clusters with number of elements less than 5.</p>
		<p><img src='GO_heatmap_random_@{ont}@{substr}_@{i}.jpg' width='1200' /></p>
		<br>
		<p><b>Figure 2.</b> Compare clustering results. Top left panel: Similarity heatmap with clusterings from different methods. Bottom left panel: Concordance between clustering methods. The concordance measures how similar two clusterings are. Right panel: The difference score, number of clusters and the block mean of different clusterings.</p>
		<p><img src='GO_heatmap_random_@{ont}@{substr}_@{i}_barplot.png' width='800' /></p>
		<br>
		<p><b>Figure 3.</b> Select proper cutoff for binary cut.</p>
		<p><img src='GO_heatmap_random_@{ont}@{substr}_@{i}_binary_cut_select_cutoff.png' width='600' /></p>
		
		</body>
		</html>
		"), con = html_file)

			m = rbind(m, sapply(lt[[2]], function(x) {
				tb = table(x)
				paste0(length(tb), " (", sum(tb >= 5), ")")
			}))
		}

		m = cbind(run = 1:100, m)
		m = cbind(m, "Details" = qq("<a href='https://simplifyEnrichment.github.io/random_@{ont}@{substr}/image/GO_heatmap_random_@{ont}@{substr}_@{1:100}.html'>view</a>", collapse = FALSE))
		
		# the main page
		writeLines(qq("
		<html>
		<head>
		<title>random_@{ont}@{substr}</title>
		<link href='main.css' rel='stylesheet'>
		</head>
		<body>
		<h1>random_@{ont}@{substr} (@{info})</h1>
		<hr>
		<p><b>Figure 1.</b>Compare clustering results. Left panel: The difference score, number of clusters and the block mean of different clusterings. Right panel: Concordance between clustering methods. The concordance measures how similar two clusterings are. The definition of the concordance score can be found <a href='/concordance.html'>here</a>.</p>
		<p><img src='random_@{ont}@{substr}_results.png' width='800' /></p>
		<p><b>Table 1.</b>Number of clusters identified by each clustering method. Numbers in the table indicate the number of clusters. The numbers inside the parentheses are the number of clusters with size >= 5. </p>
		@{knitr::kable(m, row.names = FALSE, format = 'html', escape = FALSE)}
		</body>
		</html>
		"), con = qq("random_@{ont}@{substr}.html"))
	}
}

## generate the main figures
ont = "BP"
for(method in c("", "jaccard", "kappa", "dice", "overlap")) {

	if(method == "") {
		substr = ""
	} else {
		substr = paste0("_", method)
	}

	summarize("random_@{ont}@{substr}/rds/clt_random_@{ont}@{substr}_@{i}.rds", 100, output = qq("random_@{ont}@{substr}_results.rds"), rerun = TRUE)
}

for(ont in c("CC", "MF")) {
	substr = ""
	summarize("random_@{ont}@{substr}/rds/clt_random_@{ont}@{substr}_@{i}.rds", 100, output = qq("random_@{ont}@{substr}_results.rds"), rerun = TRUE)
}


servr::httd("/omics/groups/OE0246/internal/guz/simplifyGO_test/")

