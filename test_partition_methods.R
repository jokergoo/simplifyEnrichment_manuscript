
library(GetoptLong)

setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/test_partition_methods")

library(simplifyEnrichment)
library(bsub)
bsub_opt$temp_dir = "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/bsub_temp"

######################## random GO ###################

bsub_opt$enforce = TRUE

for(i in 1:100) {
seed = round(runif(1)*2^30)	
bsub_chunk(name = qq("test_partition_methods_@{i}"), variables = c("i", "seed"), memory = 3, hour = 0.8,
{

	library(GetoptLong)

    setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/test_partition_methods")

	library(simplifyEnrichment)

	set.seed(seed)
	qqcat("random seed: @{seed}\n")

	remove_clustering_methods(all_clustering_methods())
	register_clustering_methods(
		partition_by_pam = function(mat, ...) binary_cut(mat, partition_fun = partition_by_pam),
		partition_by_kmeanspp = function(mat, ...) binary_cut(mat, partition_fun = partition_by_kmeanspp),
		partition_by_hclust = function(mat, ...) binary_cut(mat, partition_fun = partition_by_hclust)
	)

	qqcat("====== random @{i} ========\n")
	go_id = random_GO(500)
	mat = GO_similarity(go_id)

	clt = cmp_make_clusters(mat, method = all_clustering_methods())

	jpeg(qq("image/clt_test_partition_methods_@{i}.jpg"), width = 1200, height = 300)
	cmp_make_plot(mat, clt, nrow = 1, plot_type = "heatmap")
	dev.off()

	png(qq("image/clt_test_partition_methods_@{i}_barplot.png"), width = 800*2, height = 600*2, res = 72*2)
	cmp_make_plot(mat, clt)
	dev.off()

	saveRDS(list(mat, clt, seed), file = qq("rds/clt_test_partition_methods_@{i}.rds"))

	cat("done.\n")
})
}


for(i in 1:100) {
	if(!is_job_finished(qq("test_partition_methods_@{i}"))) {
		qqcat("not finished: test_partition_methods_@{i}\n")
	}
}


servr::httd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/")

remove_clustering_methods(all_clustering_methods())
register_clustering_methods(
	partition_by_pam = function(mat, ...) binary_cut(mat, partition_fun = partition_by_pam),
	partition_by_kmeanspp = function(mat, ...) binary_cut(mat, partition_fun = partition_by_kmeanspp),
	partition_by_hclust = function(mat, ...) binary_cut(mat, partition_fun = partition_by_hclust)
)

summarize("rds/clt_test_partition_methods_@{i}.rds", 100, output = "test_partition_methods_results.rds", rerun = TRUE)

m = NULL
for(i in 1:100) {
	qqcat("reading rds/clt_test_partition_methods_@{i}.rds, @{i}/100\n")
	lt = readRDS(qq("rds/clt_test_partition_methods_@{i}.rds"))

	html_file = qq("image/clt_test_partition_methods_@{i}.html")
	link = qq("<a href='clt_test_partition_methods_@{i-1}.html' title='clt_test_partition_methods_@{i-1}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='clt_test_partition_methods_@{i+1}.html' title='clt_test_partition_methods_@{i+1}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='../test_partition_methods.html'>Back</a>")
	if(i == 1) {
		link = qq("<a href='clt_test_partition_methods_100.html' title='clt_test_partition_methods_100'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='clt_test_partition_methods_@{i+1}.html' title='clt_test_partition_methods_@{i+1}'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='../test_partition_methods.html'>Back</a>")
	} else if(i == 100) {
		link = qq("<a href='clt_test_partition_methods_@{i-1}.html' title='clt_test_partition_methods_@{i-1}'>Previous</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='clt_test_partition_methods_1.html' title='clt_test_partition_methods_1'>Next</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='../test_partition_methods.html'>Back</a>")
	}
	writeLines(qq("
<html>
<head>
<title>test_partition_methods_@{i}</title>
<link href='../main.css' rel='stylesheet'>
</head>
<body>
<h2>test_partition_methods_@{i}</h2>
<hr />
<p>@{link}</p>
<p><b>Figure 1.</b> Clusterings from different methods. Green bars represent the small clusters with number of elements less than 5.</p>
<p><img src='clt_test_partition_methods_@{i}.jpg' width='1000' /></p>
<br>
<p><b>Figure 2.</b> Compare clustering results. Top left panel: Similarity heatmap with clusterings from different methods. Bottom left panel: Concordance between clustering methods. The concordance measures how similar two clusterings are. Right panel: The difference score, number of clusters and the block mean of different clusterings.</p>
<p><img src='clt_test_partition_methods_@{i}_barplot.png' width='800' /></p>

</body>
</html>
"), con = html_file)

	m = rbind(m, sapply(lt[[2]], function(x) {
		tb = table(x)
		paste0(length(tb), " (", sum(tb >= 5), ")")
	}))
}

m = cbind(run = 1:100, m)
m = cbind(m, "Details" = qq("<a href='image/clt_test_partition_methods_@{1:100}.html'>view</a>", collapse = FALSE))


# the main page
writeLines(qq("
<html>
<head>
<title>Test partitioning methods</title>
<link href='main.css' rel='stylesheet'>
</head>
<body>
<h1>Test partitioning methods</h1>
<hr>
<p>In the process of recursively splitting the similarity matrix in binary cut algorithm, 
in each iteration step, the current matrix is partitioned into two groups using PAM as default. 
Here we compare following partitioning methds: k-means++, PAM and hierarchical clustering with 
'ward.D2' method, on 500 random GO lists.</p>

<p><b>Figure 1.</b>Compare clustering results. Left panel: The difference score, number of clusters and the block mean of different clusterings. Right panel: Concordance between clustering methods. The concordance measures how similar two clusterings are. The definition of the concordance score can be found <a href='/concordance.html'>here</a>.</p>
<p><img src='test_partition_methods_results.png' width='800' /></p>
<p><b>Table 1.</b>Number of clusters identified by each clustering method. Numbers in the table indicate the number of clusters. The numbers inside the parentheses are the number of clusters with size >= 5. </p>
@{knitr::kable(m, row.names = FALSE, format = 'html', escape = FALSE)}
</body>
</html>
"), con = qq("test_partition_methods.html"))
