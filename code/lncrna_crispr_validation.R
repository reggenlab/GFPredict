################################ validation on cell cycle cluster genes using lncRNA CRISPR genes #############
library('ggplot2')
#Load the CRISPR screeen 
crispr_lncrna_cell_cycle_gm128 <- read.csv('../data/crispr_lncrna_cell_cycle_gm128.csv')
crispr_lncrna_cell_cycle_k562 <- read.csv('../data/crispr_lncrna_cell_cycle_k562.csv')

#load predicted genes of functions related to cell cycle process
cell_cycle_tpr <- read.csv("../data/cell_cycle_functions.csv")


genes_list_cell_cycle = list() #select all the genes into this list
for (i in 1:nrow(cell_cycle_tpr)) {
row_i <- cell_cycle_tpr[i,]
 if (row_i$tpr == 100) {next} #skip if tpr is 100, because 0 genes will be predicted
row_i_non <- row_i[which(row_i != '')]
genes <- row_i_non[9:length(row_i_non)]
genes_list_cell_cycle <- append(genes_list_cell_cycle, genes)
}
genes_list_cell_cycle <- unique(unlist(genes_list_cell_cycle))

#check the genes which are present in the CRISPR screen and get their scores
model_score_gm12878 <- crispr_lncrna_cell_cycle_gm128[which(crispr_lncrna_cell_cycle_gm128$Gene_symbol %in% genes_list_cell_cycle), 'Screen_score']
set.seed(1)
#Select equal number of negatives
null_score_gm12878 <- crispr_lncrna_cell_cycle_gm128[sample(1:nrow(crispr_lncrna_cell_cycle_gm128), length(model_score_gm12878)), 'Screen_score']
message("Significance of the validation")
wilcox.test(model_score_gm12878, null_score_gm12878)
# p-value = 0.01337

## making data for plotting
model_col7_1 <- cbind(model_score_gm12878, "Cell cycle cluster genes")
null_col7_1 <- cbind(null_score_gm12878, "Random genes")

positions = c("Cell cycle cluster genes", "Random genes")
plot_data_culster1 <- as.data.frame(rbind(model_col7_1, null_col7_1)); colnames(plot_data_culster1) <- c('Z_score', 'Type')
plot_data_culster1[,1] <- as.numeric(plot_data_culster1[,1])

plot_data_culster1 <- plot_data_culster1[-which(plot_data_culster1$Z_score >= 4),]

### plotting
pdf('compare_cluster_random_lncrna_gm12878.pdf', height = 4, width = 4)
ggplot2::ggplot(plot_data_culster1, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('red2', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 8),
axis.title.y = element_text(face="bold", colour="black", size = 8),
legend.title = element_text(face="bold", size = 8),
legend.text = element_text(size=8),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 10), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) +
labs(y="Score", x = "Type", title = "Scores of cell cycle cluster predicted genes against \n random genes in lncRNA CRISPR screen",
subtitle= 'GM12878') +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
dev.off()

message("Check 'compare_cluster_random_lncrna_gm12878.pdf' for plot showing 
CRISPR score of predicted genes in GM12878 lncRNA CRISPR screen")


## similarly check in different cell line

## check the predicted genes scores in CRISPR score
model_score_k562 <- crispr_lncrna_cell_cycle_k562[which(crispr_lncrna_cell_cycle_k562$Gene_symbol %in% genes_list_cell_cycle), 'Screen_score']
set.seed(1)

# Get the scores of equal number of random genes
null_score_k562 <- crispr_lncrna_cell_cycle_k562[sample(1:nrow(crispr_lncrna_cell_cycle_k562), length(model_score_k562)), 'Screen_score']
wilcox.test(model_score_k562, null_score_k562)
#p-value = 0.007253

## data for plotting
model_col7_2 <- cbind(model_score_k562, "Cell cycle cluster genes")
null_col7_2 <- cbind(null_score_k562, "Random genes")

positions = c("Cell cycle cluster genes", "Random genes")
plot_data_culster2 <- as.data.frame(rbind(model_col7_2, null_col7_2)); colnames(plot_data_culster2) <- c('Z_score', 'Type')
plot_data_culster2[,1] <- as.numeric(plot_data_culster2[,1])

plot_data_culster2 <- plot_data_culster2[-which(plot_data_culster2$Z_score >= 4),]


pdf('compare_cluster_random_lncrna_k562.pdf', height = 4, width = 4)
ggplot2::ggplot(plot_data_culster2, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('gold1', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 8),
axis.title.y = element_text(face="bold", colour="black", size = 8),
legend.title = element_text(face="bold", size = 8),
legend.text = element_text(size=8),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 10), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) +
labs(y="Score", x = "Type", title = "Scores of cell cycle cluster predicted genes against \n random genes in lncRNA CRISPR screen",
subtitle= 'K562') +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
dev.off()

message("Check 'compare_cluster_random_lncrna_k562.pdf' for plot showing 
CRISPR score of predicted genes in K562 lncRNA CRISPR screen")