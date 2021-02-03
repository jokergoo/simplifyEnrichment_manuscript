setwd("~/manuscript/simplifyEnrichment")

library(simplifyEnrichment)

set.seed(888)
go_id = random_GO(500)
mat = GO_similarity(go_id)
cl = binary_cut(mat)

pdf("figure3.pdf", width = 9, height = 5)
ht_clusters(mat, cl, column_title = "Similarity matrix of 500 random GO terms",
	exclude_words = c("process", "regulation"), use_raster = FALSE)
dev.off()
