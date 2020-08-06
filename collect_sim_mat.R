
library(cola)
library(GetoptLong)

setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test")


################# perform gene set enrichment test ##############

# data was downloaded from https://www.ebi.ac.uk/gxa/download
all_files = scan(pipe("ls dataset/*/*-analytics.tsv"), what = "character")

diff_gene_lt = list()
for(i in seq_along(all_files)) {
	cat(i, "/", length(all_files), "\n")
	tb = fread(all_files[i]); tb = as.data.frame(tb)
	if(grepl("^ENSG\\d", tb[1, 1])) {
		ind = grep(".p-value", colnames(tb))
		for(j in ind) {
			adjp = p.adjust(tb[, j], "BH"); adjp[is.na(adjp)] = 1
			g = tb[ adjp < 0.05 , 1]
			if(length(g) > 500 && length(g) < 3000) {
				nm = paste0(gsub("-analytics.tsv$", "", basename(all_files[i])), "_", gsub(".p-value", "", colnames(tb)[j]))
				diff_gene_lt[[nm]] = g
			}
		}
	}	
}

saveRDS(diff_gene_lt, file = "rds/diff_gene_lt.rds")
diff_gene_lt = readRDS("rds/diff_gene_lt.rds")

map = map_to_entrez_id("ENSEMBL")
diff_gene_lt2 = lapply(diff_gene_lt, function(x) {
	x = map[x]
	x[!is.na(x)]
})


## apply functional enrichment for each set of differnetial genes ###
library(bsub)
library(GetoptLong)
bsub_opt$temp_dir = "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/bsub_temp"
for(nm in names(diff_gene_lt2)) {
	g = diff_gene_lt2[[nm]]
	bsub_chunk(name = qq("enrichment_@{nm}"), hour = 1, memory = 2, variables = c("nm", "g"),
	{	
		library(GetoptLong)
		res = list()
		res$GO_BP = clusterProfiler::enrichGO(
	            gene = g,
	            OrgDb = "org.Hs.eg.db",
	            ont = "BP",
	            pAdjustMethod = "BH",
	            minGSSize = 10,
	            maxGSSize = 1000,
	            pvalueCutoff  = 1,
	            qvalueCutoff  = 1)
		res$KEGG = clusterProfiler::enrichKEGG(
	                    gene = g,
	                    pAdjustMethod = "BH",
	                    minGSSize = 10,
	                    maxGSSize = 1000,
	                    pvalueCutoff  = 1,
	                    qvalueCutoff  = 1)
		res$DO = DOSE::enrichDO(
	                    gene = g,
	                    ont = "DO",
	                    pAdjustMethod = "BH",
	                    minGSSize = 10,
	                    maxGSSize = 1000,
	                    pvalueCutoff  = 1,
	                    qvalueCutoff  = 1)
		res$Reactome = ReactomePA::enrichPathway(
	                    gene = g,
	                    pAdjustMethod = "BH",
	                    minGSSize = 10,
	                    maxGSSize = 1000,
	                    pvalueCutoff  = 1,
	                    qvalueCutoff  = 1)
		# MsigDB
		library(msigdbr)
		m_df = msigdbr(species = "Homo sapiens")
		set_lt = list("C2" = c("CGP"),
			          "C3" = c("MIR:MIR_Legacy", "TFT:TFT_Legacy", "TFT:GTRD", "MIR:MIRDB"),
			          "C4" = c("CM", "CGN"),
			          "C7" = "")
		for(set in names(set_lt)) {
			for(subset in set_lt[[set]]) {
				l = m_df$gs_cat == set & m_df$gs_subcat == subset
				res[[qq("MsigDB_@{set}_@{gsub('^.*:', '', subset)}")]] = clusterProfiler::enricher(
					            gene = g, 
					            TERM2GENE = as.data.frame(m_df[l, c("gs_id", "entrez_gene")]),
					            pAdjustMethod = "BH",
			                    minGSSize = 10,
			                    maxGSSize = 1000,
			                    pvalueCutoff  = 1,
			                    qvalueCutoff  = 1)
			}
		}
		
		saveRDS(res, file = qq("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/rds/enrichment/@{nm}_enrichment.rds"))
	})
}


###################### generate similarity matrix ######################

all_nm = names(diff_gene_lt)

### number of sig terms
n_lt = list()
for(i in seq_along(all_nm)) {
	nm = all_nm[i]
	
	qqcat("reading @{i}/@{length(all_nm)}.\n")
	res = readRDS(qq("rds/enrichment/@{nm}_enrichment.rds"))

	lt = list()
	l = p.adjust(res$GO_BP@result$pvalue, "BH") < 0.05
	n_lt$GO_BP = c(n_lt$GO_BP, sum(l))

	l = p.adjust(res$DO@result$pvalue, "BH") < 0.05
	n_lt$DO = c(n_lt$DO, sum(l))

	l = p.adjust(res$KEGG@result$pvalue, "BH") < 0.05
	n_lt$KEGG = c(n_lt$KEGG, sum(l))

	l = p.adjust(res$Reactome@result$pvalue, "BH") < 0.05
	n_lt$Reactome = c(n_lt$Reactome, sum(l))

	for(sm in c("C2_CGP", "C3_MIR_Legacy", "C3_TFT_Legacy", "C3_GTRD", "C3_MIRDB", "C4_CM", "C4_CGN", "C7_")) {
		l = p.adjust(res[[qq("MsigDB_@{sm}")]]@result$pvalue, "BH") < 0.05
		n_lt[[qq("MsigDB_@{sm}")]] = c(n_lt[[qq("MsigDB_@{sm}")]], sum(l))
	}
}

library(bsub)
bsub_opt$temp_dir = "/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test/bsub_temp"


for(i in seq_along(all_nm)) {
	nm = all_nm[i]
	bsub_chunk(name = qq("generate_sim_mat_@{i}"), hour = 1, memory = 1, variables = c("nm", "all_nm", "i"),
	{ 
	setwd("/icgc/dkfzlsdf/analysis/B080/guz/simplifyGO_test")

	library(simplifyEnrichment)
	library(DOSE)
	library(GetoptLong)

	qqcat("reading @{i}/@{length(all_nm)}.\n")
	res = readRDS(qq("rds/enrichment/@{nm}_enrichment.rds"))

	lt = list()
	l = p.adjust(res$GO_BP@result$pvalue, "BH") < 0.05; qqcat("  - [@{i}/@{length(all_nm)}] @{sum(l)}/@{length(l)} GO terms.\n")
	if(sum(l) >= 100 && sum(l) <= 1000) {
		lt$GO_BP_sim = GO_similarity(res$GO_BP@result$ID[l])
		lt$GO_BP_sim_kappa = term_similarity_from_enrichResult(subset_enrichResult(res$GO_BP, l), method = "kappa")
		lt$GO_BP_sim_jaccard = term_similarity_from_enrichResult(subset_enrichResult(res$GO_BP, l), method = "jaccard")
		lt$GO_BP_sim_dice = term_similarity_from_enrichResult(subset_enrichResult(res$GO_BP, l), method = "dice")
		lt$GO_BP_sim_overlap = term_similarity_from_enrichResult(subset_enrichResult(res$GO_BP, l), method = "overlap")
	}

	l = p.adjust(res$DO@result$pvalue, "BH") < 0.05; qqcat("  - [@{i}/@{length(all_nm)}] @{sum(l)}/@{length(l)} DO terms.\n")
	if(sum(l) >= 50 && sum(l) <= 1000) {
		lt$DO_sim = DO_similarity(res$DO@result$ID[l])
		lt$DO_sim_kappa = term_similarity_from_enrichResult(subset_enrichResult(res$DO, l), method = "kappa")
		lt$DO_sim_jaccard = term_similarity_from_enrichResult(subset_enrichResult(res$DO, l), method = "jaccard")
		lt$DO_sim_dice = term_similarity_from_enrichResult(subset_enrichResult(res$DO, l), method = "dice")
		lt$DO_sim_overlap = term_similarity_from_enrichResult(subset_enrichResult(res$DO, l), method = "overlap")
	}
	
	l = p.adjust(res$KEGG@result$pvalue, "BH") < 0.05; qqcat("  - [@{i}/@{length(all_nm)}] @{sum(l)}/@{length(l)} KEGG terms.\n")
	if(sum(l) >= 25 && sum(l) <= 1000) {
		lt$KEGG_sim_kappa = term_similarity_from_enrichResult(subset_enrichResult(res$KEGG, l), method = "kappa")
		lt$KEGG_sim_jaccard = term_similarity_from_enrichResult(subset_enrichResult(res$KEGG, l), method = "jaccard")
		lt$KEGG_sim_dice = term_similarity_from_enrichResult(subset_enrichResult(res$KEGG, l), method = "dice")
		lt$KEGG_sim_overlap = term_similarity_from_enrichResult(subset_enrichResult(res$KEGG, l), method = "overlap")
	}
	
	l = p.adjust(res$Reactome@result$pvalue, "BH") < 0.05; qqcat("  - [@{i}/@{length(all_nm)}] @{sum(l)}/@{length(l)} Reactome terms.\n")
	if(sum(l) >= 25 && sum(l) <= 1000) {
		lt$Reactome_sim_kappa = term_similarity_from_enrichResult(subset_enrichResult(res$Reactome, l), method = "kappa")
		lt$Reactome_sim_jaccard = term_similarity_from_enrichResult(subset_enrichResult(res$Reactome, l), method = "jaccard")
		lt$Reactome_sim_dice = term_similarity_from_enrichResult(subset_enrichResult(res$Reactome, l), method = "dice")
		lt$Reactome_sim_overlap = term_similarity_from_enrichResult(subset_enrichResult(res$Reactome, l), method = "overlap")
	}

	for(sm in c("C2_CGP", "C3_MIR_Legacy", "C3_TFT_Legacy", "C3_GTRD", "C3_MIRDB", "C4_CM", "C4_CGN", "C7_")) {
		l = p.adjust(res[[qq("MsigDB_@{sm}")]]@result$pvalue, "BH") < 0.05; qqcat("  - [@{i}/@{length(all_nm)}] @{sum(l)}/@{length(l)} MsigDB_@{sm} terms.\n")
		if(sum(l) >= 25 && sum(l) <= 1000) {
			lt[[qq("MsigDB_@{sm}_sim_kappa")]] = term_similarity_from_enrichResult(subset_enrichResult(res[[qq("MsigDB_@{sm}")]], l), method = "kappa")
			lt[[qq("MsigDB_@{sm}_sim_jaccard")]] = term_similarity_from_enrichResult(subset_enrichResult(res[[qq("MsigDB_@{sm}")]], l), method = "jaccard")
			lt[[qq("MsigDB_@{sm}_sim_dice")]] = term_similarity_from_enrichResult(subset_enrichResult(res[[qq("MsigDB_@{sm}")]], l), method = "dice")
			lt[[qq("MsigDB_@{sm}_sim_overlap")]] = term_similarity_from_enrichResult(subset_enrichResult(res[[qq("MsigDB_@{sm}")]], l), method = "overlap")
		}
	}

	saveRDS(lt, qq("rds/sim/lt_sim_@{nm}.rds"))

	})
}

lt2 = list()
for(i in seq_along(all_nm)) {
	nm = all_nm[i]
	qqcat("reading rds/sim/lt_sim_@{nm}.rds, @{i}/@{length(all_nm)}\n")
	lt = readRDS(qq("rds/sim/lt_sim_@{nm}.rds"))

	for(onto in names(lt)) {
		if(!onto %in% names(lt2)) lt2[[onto]] = list()
		lt2[[onto]][[nm]] = lt[[onto]]
	}
}

saveRDS(lt2, file = "rds/lt2_sim_all.rds")

saveRDS(lt2$GO_BP_sim, file = "rds/GO_BP_sim.rds")
saveRDS(lt2$GO_BP_sim_kappa, file = "rds/GO_BP_sim_kappa.rds")
saveRDS(lt2$GO_BP_sim_jaccard, file = "rds/GO_BP_sim_jaccard.rds")
saveRDS(lt2$GO_BP_sim_dice, file = "rds/GO_BP_sim_dice.rds")
saveRDS(lt2$GO_BP_sim_overlap, file = "rds/GO_BP_sim_overlap.rds")
saveRDS(lt2$DO_sim, file = "rds/DO_sim.rds")
saveRDS(lt2$DO_sim_kappa, file = "rds/DO_sim_kappa.rds")
saveRDS(lt2$DO_sim_jaccard, file = "rds/DO_sim_jaccard.rds")
saveRDS(lt2$DO_sim_dice, file = "rds/DO_sim_dice.rds")
saveRDS(lt2$DO_sim_overlap, file = "rds/DO_sim_overlap.rds")
saveRDS(lt2$KEGG_sim_kappa, file = "rds/KEGG_sim_kappa.rds")
saveRDS(lt2$KEGG_sim_jaccard, file = "rds/KEGG_sim_jaccard.rds")
saveRDS(lt2$KEGG_sim_dice, file = "rds/KEGG_sim_dice.rds")
saveRDS(lt2$KEGG_sim_overlap, file = "rds/KEGG_sim_overlap.rds")
saveRDS(lt2$Reactome_sim_kappa, file = "rds/Reactome_sim_kappa.rds")
saveRDS(lt2$Reactome_sim_jaccard, file = "rds/Reactome_sim_jaccard.rds")
saveRDS(lt2$Reactome_sim_dice, file = "rds/Reactome_sim_dice.rds")
saveRDS(lt2$Reactome_sim_overlap, file = "rds/Reactome_sim_overlap.rds")
for(sm in c("C2_CGP", "C3_MIR_Legacy", "C3_TFT_Legacy", "C3_GTRD", "C3_MIRDB", "C4_CM", "C4_CGN", "C7_")) {
	saveRDS(lt2[[qq("MsigDB_@{sm}_sim_kappa")]], file = qq("rds/MsigDB_@{sm}_sim_kappa.rds"))
	saveRDS(lt2[[qq("MsigDB_@{sm}_sim_jaccard")]], file = qq("rds/MsigDB_@{sm}_sim_jaccard.rds"))
	saveRDS(lt2[[qq("MsigDB_@{sm}_sim_dice")]], file = qq("rds/MsigDB_@{sm}_sim_dice.rds"))
	saveRDS(lt2[[qq("MsigDB_@{sm}_sim_overlap")]], file = qq("rds/MsigDB_@{sm}_sim_overlap.rds"))
}
