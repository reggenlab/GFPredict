library(GFPredict)
library(ggplot2)
setwd("/home/omkarc/omkar/package_testing/GFPredict/code/")
##################################################### cell cycle ###################################################################
#Check PMID 29662178 to get the research paper published on this CRISPR screen

message("Loading the cell cycle CRISPR screen")
cell_cycle_screen <- read.csv('../data/cell_cycle_screen.csv') #Viability, cell-cycle & DNA-repair genes

#Fine tune the n_bootstraps if you wish
result1_rf_all <- predict_related_genes(cell_cycle_screen$HGNC[1:50], n_bootstrap = 10, ml_model = 'random.forest', feature_type = 'all_features')
#balanced-accuracy = 0.8834039 

#select the top 30 genes present in the screen
cell_cycle_screen_pval <- cell_cycle_screen[which(cell_cycle_screen$ew136.PVAL <= 2),]
pos_present1 <- which( result1_rf_all[[1]] %in% cell_cycle_screen_pval$HGNC )[1:50]
top_20_genes_pre1 <- result1_rf_all[[1]][pos_present1]
top20_genes1 <- unique(top_20_genes_pre1)[1:30]

## getting the scores of the predicted genes comaped to random genes
scores_list_model1 <- -1 * cell_cycle_screen_pval[which( cell_cycle_screen_pval$HGNC %in% top20_genes1 ), 'ew136.BCSCORE']
set.seed(1)
#select equal number of random genes
scores_list_null1 <- -1 * na.omit(cell_cycle_screen_pval[sample(1:nrow(cell_cycle_screen_pval), length(scores_list_model1)), 'ew136.BCSCORE'])

#significance test between predicted genes' scores and random genes' scores
wilcox.test(scores_list_model1, scores_list_null1)
#p-value = 0.0001448, n = 30

### data for plotting
model_col1 <- cbind(scores_list_model1, 'Predicted genes')
null_col1 <- cbind(scores_list_null1, 'Random genes')
positions = c("Predicted genes", "Random genes")
plot_data1 <- as.data.frame(rbind(model_col1, null_col1)); colnames(plot_data1) <- c('Z_score', 'Type')
plot_data1[,1] <- as.numeric(plot_data1[,1])

#plotting
#pdf('compare_model_null_cell_cycle.pdf', height = 4, width = 4)
ggplot2::ggplot(plot_data1, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('skyblue3', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 12),
axis.title.y = element_text(face="bold", colour="black", size = 12),
legend.title = element_text(face="bold", size = 10),
legend.text = element_text(size=10),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 14), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) +
labs(y="-1 * Z-score", x = "Type", title = "Z-score of predicted genes against \n random genes in cell cycle CRISPR screen") +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
ggsave("compare_model_null_cell_cycle.pdf", width = 4.5, height = 4.5, units = 'in')
#dev.off()

system("code compare_model_null_cell_cycle.pdf")

message("Check the file 'compare_model_null_cell_cycle.pdf' for results in CRISPR screen for cell cycle")

## loading lncRNA CRISPR screen for  in cell cycle 
crispr_lncrna_cell_cycle_gm128 <- read.csv('../data/crispr_lncrna_cell_cycle_gm128.csv')
crispr_lncrna_cell_cycle_k562 <- read.csv('../data/crispr_lncrna_cell_cycle_k562.csv')

## Check the predicted genes' score in the lncRNA CRISPR screen
lncrna_inter1 <- crispr_lncrna_cell_cycle_gm128[which(crispr_lncrna_cell_cycle_gm128$Gene_symbol %in% result1_rf_all[[1]]),]
set.seed(1)
#getting the scores for random genes to compare
lncrna_inter_null1 <- crispr_lncrna_cell_cycle_gm128[sample(1:length(crispr_lncrna_cell_cycle_gm128$Screen_score),  nrow(lncrna_inter1)), ]
lnc_score_gm12878 <- lncrna_inter1$Screen_score
lnc_score_gm12878_null <- lncrna_inter_null1$Screen_score

#significance test between predicted genes' scores and random genes' scores
wilcox.test(lnc_score_gm12878, lnc_score_gm12878_null)
#p-value = 0.005636, n = 9

## Check the predicted genes' score in the lncRNA CRISPR screen
lncrna_inter2 <- crispr_lncrna_cell_cycle_k562[which(crispr_lncrna_cell_cycle_k562$Gene_symbol %in% result1_rf_all[[1]]),]
set.seed(1)
#getting the scores for random genes to compare
lncrna_inter_null2 <- crispr_lncrna_cell_cycle_k562[sample(1:length(crispr_lncrna_cell_cycle_k562$Screen_score),  nrow(lncrna_inter2)), ]
lnc_score_k562 <- lncrna_inter2$Screen_score
lnc_score_k562_null <- lncrna_inter_null2$Screen_score

#significance test between predicted genes' scores and random genes' scores
wilcox.test(lnc_score_k562, lnc_score_k562_null)
#p-value = 0.03998, n = 9


#plot data
col1 <- cbind(lnc_score_gm12878, "Predicted genes, GM12878")
col1_null <- cbind(lnc_score_gm12878_null, "Random genes, GM12878")
plotdata1_1 <- as.data.frame(rbind(col1, col1_null))
colnames(plotdata1_1) <- c('Z_score', 'Type'); plotdata1_1[,1] <- as.numeric(plotdata1_1[,1])
positions1 = c("Predicted genes, GM12878", "Random genes, GM12878")

plotdata1_1 <- plotdata1_1[-which(plotdata1_1$Z_score >= 2),]


#plotting

ggplot2::ggplot(plotdata1_1, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('mediumpurple2', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 12),
axis.title.y = element_text(face="bold", colour="black", size = 12),
legend.title = element_text(face="bold", size = 10),
legend.text = element_text(size=10),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 14), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) +
labs(y="ncRNA gene score", x = "Type", title = "Scores of predicted genes against \n random genes in lncRNA-cell cycle CRISPR screen",
subtitle = 'GM12878') +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
ggsave("compare_model_null_cell_cycle_lncrna_gm12878.pdf", width = 4.5, height = 4.5, units = 'in')

system("code compare_model_null_cell_cycle_lncrna_gm12878.pdf")

## plot data
col2 <- cbind(lnc_score_k562, "Predicted, K562")
col2_null <- cbind(lnc_score_k562_null, "Random genes, K562")
plotdata1_2 <- as.data.frame(rbind(col2, col2_null))
colnames(plotdata1_2) <- c('Z_score', 'Type'); plotdata1_2[,1] <- as.numeric(plotdata1_2[,1])
positions2 = c("Predicted, K562", "Random genes, K562")

plotdata1_2 <- plotdata1_2[-which(plotdata1_2$Z_score >= 2),]


#plotting
ggplot2::ggplot(plotdata1_2, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('mediumpurple2', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 12),
axis.title.y = element_text(face="bold", colour="black", size = 12),
legend.title = element_text(face="bold", size = 10),
legend.text = element_text(size=10),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 14), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) +
labs(y="ncRNA gene score", x = "Type", title = "Scores of predicted genes against \n random genes in lncRNA-cell cycle CRISPR screen",
subtitle = 'K562') +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
ggsave("compare_model_null_cell_cycle_lncrna_k562.pdf", width = 4.5, height = 4.5, units = 'in')



system("code compare_model_null_cell_cycle_lncrna_k562.pdf")

message("Check files 'compare_model_null_cell_cycle_lncrna_gm12878', 'compare_model_null_cell_cycle_lncrna_gm12878' ")
message("For model's performance in predicting lncRNA genes after training on coding genes")




################################################ pepetide accumulation ########################################################################
message("Loading CRISPR screen for peptide accumulation")
#Check PMID 30449619 for the article related to this CRISPR screen, ew1310
peptide_accu_screen <- read.csv("../data/peptide_accu_screen.csv")

## tain the model 
result2_rf <- predict_related_genes(peptide_accu_screen$HGNC[1:50], n_bootstrap = 7, ml_model = 'random.forest', feature_type = 'all_features')
#balanced-accuracy =  0.8130740
message(paste("balanced accuracy:", result2_rf[[2]][[4]]))

#select the top 30 genes present in the screen
pos_present2 <- which( result2_rf[[1]] %in% peptide_accu_screen$HGNC )[1:50]
top_20_genes_pre2 <- result2_rf[[1]][pos_present2]
top20_genes2 <- unique(top_20_genes_pre2)[1:30]

#taking the socres of the predicted genes
scores_list_model2 <- na.omit(peptide_accu_screen[which( peptide_accu_screen$HGNC %in% top20_genes2), 4])
scores_list_model2 <- -1 * as.numeric(scores_list_model2)
set.seed(1)
#Getting the scores of random genes to compare
scores_list_null2 <-  na.omit(peptide_accu_screen[sample(1:nrow(peptide_accu_screen), length(scores_list_model2)), 4])
scores_list_null2 <- -1 * as.numeric(scores_list_null2)

#significance test between predicted genes' scores and random genes' scores
wilcox.test(scores_list_model2, scores_list_null2)
message(paste("model, median:", median(scores_list_model2),  "null model, median:", median(scores_list_null2)))

#p-value = 4.545e-06, n = 30

#data for plotting
model_col2 <- cbind(scores_list_model2, 'Predicted genes')
null_col2 <- cbind(scores_list_null2, 'Random genes')
positions = c("Predicted genes", "Random genes")
plot_data2 <- as.data.frame(rbind(model_col2, null_col2)); colnames(plot_data2) <- c('Z_score', 'Type')
plot_data2[,1] <- as.numeric(plot_data2[,1])


ggplot(plot_data2, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('skyblue3', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 12),
axis.title.y = element_text(face="bold", colour="black", size = 12),
legend.title = element_text(face="bold", size = 10),
legend.text = element_text(size=10),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 14), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) + 
labs(y="Z-score", x = "Type", title = "Z-score of predicted genes against \n random genes in peptide accumulation CRISPR screen") +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
ggsave("compare_model_null_peptide_accu.pdf", width = 4.5, height = 4.5, units = 'in')

system("code compare_model_null_peptide_accu.pdf")
 
message("Check  'compare_model_null_peptide_accu.pdf' to check the model's performance on CRISPR screen for peptide accumulation")

################################################ Resistance to chemicals ########################################################################
message("Loading data related to the Resistance to chemicals")
#Check PMID 31761535 to get the article related to this screen")
resis_chemicals_screen <- read.csv('/home/omkarc/omkar/CRISPR/all_crispr/ew1367.PubMed.gene_bcscore.csv', sep = "\t")
write.csv(resis_chemicals_screen, file = "../data/resis_chemicals_screen.csv", row.names = F)

resis_chemicals_screen <- read.csv("../data/resis_chemicals_screen.csv", )

#Run the model
result3_rf <- predict_related_genes(resis_chemicals_screen$HGNC[1:50], n_bootstrap = 9, ml_model = 'random.forest')
#balanced-accuracy = 0.8183869

#select the top 30 genes present in the screen
pos_present3 <- which( result3_rf[[1]] %in% resis_chemicals_screen$HGNC)
top_20_genes3 <- na.omit(result3_rf[[1]][pos_present3][1:30])

### Take the scores of the predicted genes
scores_list_model3 <- -1 * na.omit(resis_chemicals_screen[which( resis_chemicals_screen$HGNC %in% top_20_genes3 ), 4])
set.seed(1)
#Getting the scores of random genes to compare
scores_list_null3 <-  -1 * na.omit(resis_chemicals_screen[sample(1:nrow(resis_chemicals_screen), length(scores_list_model3)), 4])

#significance test between predicted genes' scores and random genes' scores
wilcox.test(scores_list_model3, scores_list_null3)
#p-value = 0.005082 n = 18
paste("median, model:", median(scores_list_model3), "median, null model", median(scores_list_null3))

#data for plotting
model_col3 <- cbind(scores_list_model3, 'Predicted genes')
null_col3<- cbind(scores_list_null3, 'Random genes')
positions = c("Predicted genes", "Random genes")
plot_data3 <- as.data.frame(rbind(model_col3, null_col3)); colnames(plot_data3) <- c('Z_score', 'Type')
plot_data3[,1] <- as.numeric(plot_data3[,1])

ggplot2::ggplot(plot_data3, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('skyblue3', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 12),
axis.title.y = element_text(face="bold", colour="black", size = 12),
legend.title = element_text(face="bold", size = 10),
legend.text = element_text(size=10),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 14), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) +
labs(y="-1 * Z-score", x = "Type", title = "Z-score of predicted genes against \n random genes in resistance to chemicals CRISPR screen") +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
ggsave("compare_model_null_resis_chemicals.pdf", width = 4.5, height = 4.5, units = 'in')

system("code compare_model_null_resis_chemicals.pdf")
    message("Check  'compare_model_null_resis_chemicals.pdf' to check the model's performance on CRISPR screen for resistance to chemicals")

################################################ viabiility ########################################################################
message("Loading data related to the resistance to chemicals CRISPR")
#Check PMID 31316073 to get the article related to this screen")
viability_screen <- read.csv("../data/viability_screen.csv")

#training the model
result4_rf <- predict_related_genes(viability_screen$HGNC[1:50], n_bootstrap = 6, ml_model = 'random.forest')
# balanced-accuracy = 0.8208397

#selecting the top 30 genes
pos_present4 <- which(result4_rf[[1]] %in% viability_screen$HGNC)
top_20_genes4 <- na.omit(result4_rf[[1]][pos_present4][1:30])

#getting the scores of the predicted genes
scores_list_model4 <- -1 * na.omit(viability_screen[which(viability_screen$HGNC %in% top_20_genes4), 4])
set.seed(1)
#Getting the scores of random genes to compare
scores_list_null4 <-  -1 * na.omit(viability_screen[sample(1:nrow(viability_screen), length(scores_list_model4)), 4])

#significance test between predicted genes' scores and random genes' scores
wilcox.test(scores_list_model4, scores_list_null4)
#p-value = 0.00009532, n = 30
paste("median, model:", median(scores_list_model4), "median, null model", median(scores_list_null4))

#data for plotting
model_col4 <- cbind(scores_list_model4, 'Predicted genes')
null_col4 <- cbind(scores_list_null4, 'Random genes')
positions = c("Predicted genes", "Random genes")
plot_data4 <- as.data.frame(rbind(model_col4, null_col4)); colnames(plot_data4) <- c('Z_score', 'Type')
plot_data4[,1] <- as.numeric(plot_data4[,1])

#plotting

ggplot2::ggplot(plot_data4, aes(x = Type, y = Z_score, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15, width = 0.3) +
scale_fill_manual(values=c('skyblue3', 'grey55')) +
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 8),
axis.title.x = element_text(face="bold", colour="black", size = 12),
axis.title.y = element_text(face="bold", colour="black", size = 12),
legend.title = element_text(face="bold", size = 10),
legend.text = element_text(size=10),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 14), axis.text.x=element_blank(),
axis.ticks.x=element_blank(), legend.position = 'bottom', 
legend.direction = 'horizontal') +
scale_color_discrete(labels = positions) +
labs(y="-1 * Z-score", x = "Type", title = "Z-score of predicted genes against \n random genes in viability CRISPR screen") +
guides(fill=guide_legend(nrow=1,byrow=TRUE)) 
ggsave("compare_model_null_viability.pdf", width = 4.5, height = 4.5, units = 'in')

system("code compare_model_null_viability.pdf")
message("Check  'compare_model_null_viability.pdf' to check the model's performance on CRISPR screen for viability")
