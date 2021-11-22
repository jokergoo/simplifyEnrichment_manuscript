setwd("~/manuscript/simplifyEnrichment/figures")

lt_foo = readRDS("random_BP_results.rds");
diff_score_df = lt_foo$diff_score_df
cluster_number_df = lt_foo$cluster_number_df
cluster_number2_df = lt_foo$cluster_number2_df
block_mean_df = lt_foo$block_mean_df
	

library(ggplot2)
library(cowplot)

p1 = ggplot(diff_score_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Difference score") + ggtitle("A) On random GO lists") +
	theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) 
cluster_number_df$type = "All sizes"
cluster_number2_df$type = "Size >= 5"
p2 = ggplot(rbind(cluster_number_df, cluster_number2_df), aes(x = method, y = value, col = type)) +
	geom_boxplot() + ylab("Numbers of clusters") + labs(col = "") + ggtitle("B) On random GO lists") +
	theme(axis.title.x = element_blank() ,axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(legend.position = c(0.99, 0.99), legend.justification = c("right", "top"), legend.direction="horizontal")
p3 = ggplot(block_mean_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Block mean") + ggtitle("C) On random GO lists") +
	theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))


lt_foo = readRDS("EBI_Expression_Atlas_GO_BP_results.rds"); 
diff_score_df = lt_foo$diff_score_df
cluster_number_df = lt_foo$cluster_number_df
cluster_number2_df = lt_foo$cluster_number2_df
block_mean_df = lt_foo$block_mean_df
	

library(ggplot2)
library(cowplot)

p4 = ggplot(diff_score_df, aes(x = method, y = value)) +
	geom_boxplot() + ylab("Difference score") + ggtitle("D) On Expression Atlas datasets") +
	theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
cluster_number_df$type = "All sizes"
cluster_number2_df$type = "Size >= 5"
p5 = ggplot(rbind(cluster_number_df, cluster_number2_df), aes(x = method, y = value, col = type)) +
	geom_boxplot() + ylab("Numbers of clusters") + labs(col = "") + ggtitle("E) On Expression Atlas datasets") +
	theme(axis.title.x = element_blank() ,axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(legend.position = c(0.99, 0.99), legend.justification = c("right", "top"), legend.direction="horizontal")
p6 = ggplot(block_mean_df, aes(x = method, y = value)) + ggtitle("F) On Expression Atlas datasets") +
	geom_boxplot() + ylab("Block mean") +
	theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

pdf("figure5.pdf", width = 12, height = 8)
print(plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2))
dev.off()



