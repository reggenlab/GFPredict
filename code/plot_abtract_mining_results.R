library(ggplot2)
    
term_freq <- read.csv('../data/term_freq.csv') # reading pre-run data for predicted genes containing number of co-occrrences

null_freq <- read.csv('../data/null_freq.csv') # reading pre-run data for random genes containing number of co-occrrences


get_plot_values <- function (data) {

data_rm_na <- data[, colSums(is.na(data)) != nrow(data)]
count_list <- list()
for (row in 1:nrow(data_rm_na)) {
    row_i <- data_rm_na[row,]
    count_i <- length(which(row_i[5:length(row_i)] != 'NA'))
    count_list <- append(count_list, count_i) 
}
count_list <- unlist(count_list)
}

count_for_predicted_genes_function_term <- get_plot_values(term_freq)
count_for_random_genes_function_term <- get_plot_values(null_freq)


## Add zeros to the random gene list so that it will be of equal length to predicted genes list
zeros_to_add <- length(count_for_predicted_genes_function_term) - length(count_for_random_genes_function_term)
zero_list <- list(rep(0, zeros_to_add) )
count_for_random_genes_function_term_zeros <- c(count_for_random_genes_function_term, zero_list[[1]]   )

#creating plot data
counts_predicted <- cbind(count_for_predicted_genes_function_term, "Predicted genes")
counts_random <- cbind(count_for_random_genes_function_term_zeros, "Random genes")
plot_data <- as.data.frame(rbind(counts_predicted, counts_random))
colnames(plot_data) <- c('count', 'Type')
plot_data[,1] <- as.numeric(plot_data[,1])

pdf('box_plot_cooccurrence.pdf', height = 4, width = 5)
ggplot(plot_data, aes(x = Type, y = count, fill = Type)) +
geom_boxplot(lwd = 0.2, outlier.size = 0.15) +
scale_fill_manual(values=c('dodgerblue3', 'darkslategray')) + 
theme(plot.title = element_text(hjust = 0.5, face="bold", colour="black", size = 10), 
axis.title.x = element_text(face="bold", colour="black", size = 10),
axis.title.y = element_text(face="bold", colour="black", size = 10),
legend.title = element_text(face="bold", size = 10),
panel.background = element_rect(fill='white', colour='white'),
axis.text.y = element_text(hjust = +1, size = 12),
axis.text.x = element_text(size = 12)) + 
labs(y="Function term and gene term \n co-occurrence count", x = "Model type", title = 'Co-occurrence of predicted gene term with \n corresponding function term 
in PubMed abstracts') 
dev.off()

message("Check the pdf file 'box_plot_cooccurrence.pdf' to see the boxplot comparing the co-occurrence counts between model (predicted genes) and null modell (random genes).")
