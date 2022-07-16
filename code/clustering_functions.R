library(Rtsne)
suppressPackageStartupMessages(library(factoextra))


message('Loading data ...')
tf_mat_tpr <- read.csv('../data/tf_mat_tpr_clustering.csv') 
#loading top predictors's correlation with the predictions of random forest model for multiple functions whose confidence score or tpr is above 60
tf_cor_tpr <- read.csv('../data/tf_cor_tpr_clustering.csv')

#This below function calculates the score between each function based on shared top predictors
#
setoverlaps <- function(tf_mat_tpr, tf_cor_tpr) { 
nump = ncol(tf_mat_tpr) ;
jaccards = matrix( 0 , nump , nump) ;
for (i in 1:nump) {
for( j in i:nump) {
    #print(i)
    #print(j)
gs1 = tf_mat_tpr[,i] ; gs2 = tf_mat_tpr[,j] ;
pos = which(gs1 != "") ; pos1 = which(gs2 != "") ;
gs1 = gs1[pos] ; gs2 = gs2[pos1] ;
score = 0
for (t1 in gs1) {
    for (t2 in gs2) {
    if (t1 == t2) {

    if (tf_cor_tpr[which(gs1 %in% t1),i] *  tf_cor_tpr[which(gs2 %in% t2),j] > 0) {score = score + 1} # +1 score if two top predictors/TFs have same direction
    if (tf_cor_tpr[which(gs1 %in% t1),i] *  tf_cor_tpr[which(gs2 %in% t2),j] < 0) {score = score - 1} # -1 score if two top predictors/TFs have oppositte direction


    }
    }
        } 
jaccards[i,j] = score
} }

return(jaccards) ;
}
message('building similarity matrix...')
jaccard_index <- setoverlaps(as.matrix(tf_mat_tpr), as.matrix(tf_cor_tpr)) #run the above function on the data
colnames(jaccard_index) <- colnames(tf_mat_tpr) 
rownames(jaccard_index) <- colnames(tf_mat_tpr)

diag(jaccard_index) <- 0
jaccard_index[jaccard_index <= 2] <- 0
pos_zero <- which(rowSums(jaccard_index) == 0)
jaccard_zero <- jaccard_index[-pos_zero, -pos_zero]

set.seed(3) 
tsne_out <- Rtsne((10 - jaccard_zero), is_distance=T, perplexity = 50) #Run Rtsne
tsne_ordinates <-  as.data.frame(tsne_out$Y)
Dbscan_cl1 <- dbscan::dbscan(tsne_out$Y, eps = 0.16, minPts = 4) #Get clusters using DBScan
Dbscan_cl2 <- fpc::dbscan(tsne_out$Y, eps = 0.16, MinPts = 4) #Confirm if the same clusters rendered from other mehtod
check_if_true <- all(Dbscan_cl1$cluster == Dbscan_cl2$cluster)
tsne_ordi <- tsne_ordinates
tsne_ordi$cluster = Dbscan_cl1$cluster
rownames(tsne_ordi) <- rownames(jaccard_zero)
message("cluster with cell_cycle related functions")
rownames(tsne_ordi)[tsne_ordi$cluster == 47] #clusters related to cell cycle processes as mentioned in the research paper


pdf('cluster_of_functions.pdf', width=20, height=20) #This plots the image containing clusters
factoextra::fviz_cluster(Dbscan_cl1, tsne_ordinates, repel = FALSE, pointsize = 2, main = "Cluster of functions based on shared top predictors",
xlab = 'PC-1', ylab = 'PC-2', labelsize = 0, shape = "", outlier.pointsize = 1.2, geom = "point") +
geom_point() + 
theme(legend.position = "none",panel.background = element_rect(fill='white', colour='white'), axis.title=element_text(size=20,face="bold"),
axis.text = element_text(size = 16), plot.title = element_text(size=22, face = "bold", hjust = 0.5)) 
dev.off()

message("open cluster_of_functions.pdf file to view the clusters")