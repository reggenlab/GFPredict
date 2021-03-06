options(warn=-1);
options(scipen = 999);

# list.of.packages <- c("ggplot2", "Rcpp", "randomForest", "mltest", 'ROCR', "glmnet", "e1071", "xgboost")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)


#' The function makes gene prediction based on similarity in shared transcription factors 
#' and epigenomic marks present across promoter sites of the genes
#' @param genes Object of class list, containing functionally related genes
#' @param ml_model Object of class string, Machine learning model to select. Options:
#' 'xg.boost' , 'svm', 'random.forest' (default), 'linear.regression', 'logistic.regression'
#' @param n_bootstrap Number of bootstraps of the negative samples (n = 3 by default)
#' @param feature_type Object of class string. Type of features to train the model. Options: 'all_features' (default), 'transcription_factors'
#' @examples
#' cell_cycle_related_genes <- c("RPL23A", "IPO4", "TARS1", "COG4", "TEX10", "TRMT112", "TXNL4A", "CLP1", 
#' "HSPA5", "ANAPC4", "RNF168", "ATP5F1D", "RUVBL2", "NMT1", "PNKP", "SF3B3", "FDPS", "FARSB", "HARS1",
#' "RPL8", "CCT3", "AQR", "MCL1", "CENPM") 
#'  Reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5642621/
#' 
#' result <- predict_related_genes(genes = cell_cycle_related_genes)
#' 
#' result <- predict_related_genes(genes = cell_cycle_related_genes, n_bootstrap = 5)
#' 
#' result <- predict_related_genes(genes = cell_cycle_related_genes, ml_model = 'svm', n_bootstrap = 10, feature_type = 'transcription_factors')
#' 
#' result[[1]] # contains list of predicted genes 
#' 
#' result[[2]] # contains table of model performance metrics 
#' 
#' result[[3]] # contains table of top predictors (if ml_model = random.forest)
#' 
#' @return List containing predicted gene, top predictors (if ml_model = 'random forest') and ML model performance metrics
#' @export
predict_related_genes <- function(genes, ml_model, n_bootstrap, feature_type){

    if (length(unique(genes)) <= 20) {stop("Number of genes must at least 20")}
    if (missing(feature_type)) {feature_type = 'all_features' ; message("Default feature_type: 'all_features' selected")}
    if (missing(ml_model)) {ml_model = 'random.forest'; message("Default 'random.forest' model is selected")}
    if (missing(n_bootstrap)) {n_bootstrap = 3; message("Default n_bootstrap iterations is set to 3")}
    if (n_bootstrap < 2) {stop("Bootstrap n must be more than 1")}
    if (ml_model == 'xg.boost' && n_bootstrap < 15) {
        message("Recommended n_bootstrap is at least 15 for XGBoost model") }
    if (ml_model == 'xg.boost' && length(genes) < 60) {
        message("For XGBoost model recommended input genes must be more than 60")}
    if (ml_model == 'svm' ) {
        message("This will take few minutes");
        if( n_bootstrap < 5) { message("recommended n_bootstrap is at least 5 for SVM model")}
        }


present0 <- grep('peakscores_all', ls(envir=.GlobalEnv))
present1 <- grep('peakscores_tf', ls(envir=.GlobalEnv))

if (length(present0) == 0 || length(present1) == 0) {
    message('Loading data...')

meta_peak_name_tissue_path <<- system.file("extdata", "meta_name_tissue_peakscores.csv", package = "GFPredict")
meta_peak_name_tissue <<- as.matrix(data.table::fread(meta_peak_name_tissue_path))

meta_peak_path <<- system.file("extdata", "meta_data_peakscores.csv", package = "GFPredict")
meta_peak <<- as.matrix(data.table::fread(meta_peak_path, header = T))

unionPeaks_path <<- system.file("extdata", "unionpeak2.0", package = "GFPredict")
unionPeaks <<- as.matrix(data.table::fread(unionPeaks_path, header=F));

peakscores_all <<- matrix(NA, 89747, 0)
for (p_file in 1:42) {
file_name = paste('peakall', p_file, '.rds', sep = "") 
peakscore_path <- system.file("extdata", file_name, package = "GFPredict")
peakscore_n <- readRDS(peakscore_path)
peakscores_all <<- cbind(peakscores_all, peakscore_n)
}

    message('...')

peakscores_tf <<- matrix(NA, 89747, 0)
for (p_file in 1:18) {
file_name = paste('peaktf', p_file, '.rds', sep = "") 
peakscore_path <- system.file("extdata", file_name, package = "GFPredict")
peakscore_n <- readRDS(peakscore_path)
peakscores_tf <<- cbind(peakscores_tf, peakscore_n)
}
}

if (feature_type == 'all_features'){
    peakscores <- peakscores_all}

if (feature_type == 'transcription_factors'){
    peakscores <- peakscores_tf }


updim = dim(unionPeaks) ;
pos = which( unionPeaks[,4] %in% genes ) ;
posScores = peakscores[pos,] ;
posdim = dim(posScores) ;

aTPR = matrix(NA , n_bootstrap , 1000) ;
aFPR = matrix(NA , n_bootstrap , 1000) ;
aresults = matrix(NA , n_bootstrap , 1000) ;
alltfs <- matrix(NA, 20, n_bootstrap) ; 
allpred = matrix(0, updim[1], n_bootstrap); 
alltfs_score <- matrix(NA, 20, n_bootstrap) ;
alltf_corr <- matrix(NA, 20, n_bootstrap); 
impvar_mat <- matrix(NA, ncol(peakscores), n_bootstrap); 
rownames(impvar_mat) <- colnames(peakscores)
sen_spe_mat_i <- matrix(NA, n_bootstrap, 9) ;
colnames(sen_spe_mat_i) <- c("Iterations", "cutoff", "accuracy", "balanced-accuracy", "recall/sensitivity", "specificity", "F1", "MCC", "error-rate") ; 
 

calculate_results <- function(result, you_test) {  tryCatch({

predicate <- ROCR::prediction(result, yout_test)
roc.perf = ROCR::performance(predicate, measure = "tpr", x.measure = "fpr")
perf <- ROCR::performance(predicate,"tpr","fpr")
opt.cut = function(perf, predicate){
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, predicate@cutoffs)
}

cutoff <- c(iter,opt.cut(roc.perf, predicate ))
predicted_result <- ifelse(result>cutoff[4], 1, 0)
x_report <- mltest::ml_test(predicted_result, yout_test) ; 
report <- c(n_bootstrap, cutoff[4], x_report$accuracy, x_report$balanced.accuracy[2], x_report$recall[2], x_report$specificity[2], x_report$F1[2], x_report$MCC[2], x_report$error.rate)
report <- as.matrix(report) ; report <- t(report)
sen_spe_cut <- rbind(sen_spe_cut, report)
return(sen_spe_cut)
},
 error=function(e){cat("ERROR :",conditionMessage(e), "\n") })

}


for (iter in 1:n_bootstrap) {
sen_spe_cut <- matrix(NA, 0, 9) ; 
colnames(sen_spe_cut) <- c("Iterations", "cutoff", "accuracy", "balanced-accuracy", "recall/sensitivity", "specificity", "F1", "MCC", "error-rate") ; 
set.seed(iter)
negRows = sample(1:updim[1], length(pos), 1) ;
negScores = peakscores[as.numeric(negRows),] ;
negScores = peakscores[negRows,] ;
ascores = rbind(posScores , negScores) ;
ascores = as.matrix(ascores);
yout = matrix(0, (posdim[1] * 2), 1) ;
yout[1:posdim[1]] = 1 ;
samp_yout <- floor(0.75 * nrow(yout))
train_yout <- sample(seq_len(nrow(yout)), size = samp_yout); train_yout <- as.matrix(train_yout) #taking the sample of indices for the train
test_yout =  sample((1:nrow(yout))[-train_yout]); #taking rest of the it for test
yout_test <- yout[test_yout]; yout_test <- as.numeric(yout_test) # this is taking target variables
ascores_test <- ascores[test_yout,] ;
yout_training <- yout[train_yout]; yout_training <- as.numeric(yout_training)
ascores_training <- ascores[train_yout,]


################################################################################################################################################
###########################     linera regression    ##################################       
################################################################################################################################################

if(ml_model == 'linear.regression') {
	
yout_test <- data.matrix(yout_test)
yout_training <- data.matrix(yout_training)
lm_train2 <- glmnet::cv.glmnet(ascores_training , yout_training, alpha = 1)
glm_train <- glmnet::glmnet(ascores_training, yout_training, alpha = 1, lambda = lm_train2$lambda.min)

newx = data.matrix(ascores_test)
result = predict( newx= ascores_test, glm_train, s='lambda.min') ;

finalpred = predict( glm_train , data.matrix(peakscores)) ;
sidx = order(result, decreasing=TRUE) ;labels <- yout_test[sidx]
TPR=cumsum(labels)/sum(labels); 
FPR= cumsum(!labels)/sum(!labels) ;
aTPR[ iter,1:length(TPR) ] = TPR ;
aFPR[ iter,1:length(FPR) ] = FPR ;
aresults[ iter,1:length(result) ] = result[iter] ;
allpred[,iter] = finalpred ;
sen_spe_i <- calculate_results(result, yout_test)
if (length(sen_spe_i) == 0) {stop("Please try a different model (recommended: 'random.forest' or 'svm')")}
sen_spe_mat_i[iter,] <-  sen_spe_i
}

################################################################################################################################################
###########################    Logistic regression    ##################################       
################################################################################################################################################
if (ml_model == 'logistic.regression') {

yout_test <- data.matrix(yout_test)
yout_test <- as.numeric(yout_test) # nolint
yout_training <- data.matrix(yout_training)
ascores_training <- data.matrix(ascores_training)
lm_train1 <- glmnet::cv.glmnet(ascores_training, yout_training, alpha = 0, family = "binomial")
glm_train <- glmnet::glmnet(ascores_training, yout_training, alpha = 0, lambda = lm_train1$lambda.min, family = "binomial")

newx <- data.matrix(ascores_test)
result <- predict(newx = ascores_test, glm_train, s = "lambda.min")
finalpred <- predict(glm_train, data.matrix(peakscores))
sidx <- order(result, decreasing = TRUE)
labels <- yout_test[sidx]
TPR <- cumsum(labels) / sum(labels)
FPR <- cumsum(!labels) / sum(!labels)
aTPR[iter, 1:length(TPR)] <- TPR
aFPR[iter, 1:length(FPR)] <- FPR
aresults[iter, 1:length(result)] <- result[sidx]
allpred[, iter] <- finalpred
sen_spe_i <- calculate_results(result, yout_test)
if (length(sen_spe_i) == 0) {stop("Please try a different model (recommended: 'random.forest' or 'svm')")}
sen_spe_mat_i[iter,] <-  sen_spe_i

}


################################################################################################################################################
###########################     random forest model     ##################################       
################################################################################################################################################


if (ml_model == 'random.forest') {

lm_train3 <- randomForest::randomForest(yout_training ~., ntree = 50 , data=data.frame(ascores_training)) ;
#lm_train3 <- cforest(yout_training ~., data=data.frame(ascores_training),  ntree = 500L)
predtrain3 = predict(lm_train3, data.frame(ascores_training)) ;

result = predict( lm_train3 , data.frame(ascores_test)) ;

finalpred = predict( lm_train3 , data.frame(peakscores)) ;
sidx = order(result, decreasing=TRUE) ;
labels <- yout_test[sidx]; labels <- as.numeric(labels)
TPR=cumsum(labels)/sum(labels); FPR= cumsum(!labels)/sum(!labels) ;
aresults[ iter,1:length(result) ] = result[sidx] ;
aTPR[ iter,1:length(TPR) ] = TPR ;          
aFPR[ iter,1:length(FPR) ] = FPR ;
allpred[,iter ] = finalpred ;
sen_spe_i <- calculate_results(result, yout_test)
if (length(sen_spe_i) == 0) {stop("Please try adding more number of input genes")}
sen_spe_mat_i[iter,] <-  sen_spe_i

imp <- as.matrix(randomForest::importance(lm_train3, scale = TRUE ))
impvar_mat[rownames(imp),iter] <- as.numeric(imp[,1])

}



################################################################################################################################################
###########################     SVM    ##################################       
################################################################################################################################################

if (ml_model == 'svm') {

lm_train4 <- e1071::svm(yout_training ~ . , data= data.frame(ascores_training), verbose = F)
predtrain4 = predict( lm_train4 , data.frame(ascores_training)) ;

result = predict( lm_train4 , data.frame(ascores_test)) ;
finalpred = predict( lm_train4, data.frame(peakscores)) ;
sidx = order(result, decreasing=TRUE) ;
labels <- yout_test[sidx]
TPR=cumsum(labels)/sum(labels); FPR= cumsum(!labels)/sum(!labels) ;
aTPR[ iter,1:length(TPR) ] = TPR ;
aFPR[ iter,1:length(FPR) ] = FPR ;
aresults[ iter,1:length(result) ] = result[iter] ;
allpred[,iter ] = finalpred ;
sen_spe_i <- calculate_results(result, yout_test)
if (length(sen_spe_i) == 0) {stop("Please try adding more number of input genes or use ml_model = 'random.forest'")}
sen_spe_mat_i[iter,] <-  sen_spe_i


}

################################################################################################################################################
###########################    XGBoost   ##################################       
################################################################################################################################################

if (ml_model == 'xg.boost') {

samp_yout<- floor(0.75 * nrow(yout))
train_yout<-sample(seq_len(nrow(yout)), size = samp_yout); train_yout <- as.matrix(train_yout)
test_yout =  (1:nrow(yout))[-train_yout];

yout_test <- yout[test_yout]
ascores_test <- ascores[test_yout,] ;

train_id <- sample(1:length(train_yout), size = floor(0.8 * length(train_yout)), replace=FALSE) 

yout_training_id <- (train_yout[train_id, ]);   
yout_training <- yout[yout_training_id ]
ascores_training <- ascores[yout_training_id,];

yout_valid_id <- train_yout[-train_id,]; 
yout_validation <- yout[yout_valid_id,]
ascores_validation <- ascores[yout_valid_id,]

yout_training <- as.matrix(yout_training);
ascores_training <- as.matrix(ascores_training)
ascores_test <- as.matrix(ascores_test)
yout_test <- as.matrix(yout_test)
ascores_validation <- as.matrix(ascores_validation)
yout_validation <- as.matrix(yout_validation)

dtrain <- xgboost::xgb.DMatrix(ascores_training, label = (yout_training))
dvalid <- xgboost::xgb.DMatrix(data = ascores_validation, label = yout_validation)
dtest <- xgboost::xgb.DMatrix(ascores_test, label = yout_test)

lowest_error_list = list()
parameters_list = list()


for (iter_xg in 1:10){
  param <- list(booster = "gbtree",
                objective = "binary:logistic",
                max_depth = sample(2:8, 1),
                eta = runif(1, .01, .3),
                subsample = runif(1, .7, 1),
                colsample_bytree = runif(1, .8, 1),
                min_child_weight = sample(0:10, 1)
  )
  parameters <- as.data.frame(param)
  parameters_list[[iter_xg]] <- parameters
}

parameters_df = do.call(rbind, parameters_list)

for (row in 1:nrow(parameters_df)){
  set.seed(20)
  mdcv <- xgboost::xgb.train(data=dtrain,
                    booster = "gbtree",
                    objective = "binary:logistic",
                    max_depth = parameters_df$max_depth[row],
                    eta = parameters_df$eta[row],
                    subsample = parameters_df$subsample[row],
                    colsample_bytree = parameters_df$colsample_bytree[row],
                    min_child_weight = parameters_df$min_child_weight[row],
                    nrounds= 300,
                    n_jobs = 10,
                    eval_metric = "error",
                    early_stopping_rounds= 30,
                    print_every_n = 100,
                    verbose = 0,
                    watchlist = list(train= dtrain, val= dvalid)
  )
  lowest_error <- as.data.frame(1 - min(mdcv$evaluation_log$val_error))
  lowest_error_list[[row]] <- lowest_error
}

# Create object that contains all accuracy's
lowest_error_df = do.call(rbind, lowest_error_list)

# Bind columns of accuracy values and random hyperparameter values
randomsearch = cbind(lowest_error_df, parameters_df)
params <- list(booster = "gbtree", 
               objective = "binary:logistic",
               max_depth = randomsearch[1,]$max_depth,
               eta = randomsearch[1,]$eta,
               subsample = randomsearch[1,]$subsample,
               colsample_bytree = randomsearch[1,]$colsample_bytree,
               min_child_weight = randomsearch[1,]$min_child_weight)
lm_train5 <- xgboost::xgb.train(params = params,
                       data = dtrain,
                       nrounds =1000,
                       print_every_n = 10,
		       n_jobs = 10,
                       eval_metric = "auc",
                       eval_metric = "error",
                       early_stopping_rounds = 30,
                       verbose = 0,
                       watchlist = list(train= dtrain, val= dvalid))

result <- predict(lm_train5, as.matrix(ascores_test));
if (length(unique(result)) == 1) {stop("Please try a different model (recommended: 'random.forest' or 'svm') ")}
finalpred <- predict(lm_train5, as.matrix(peakscores));
sidx <- order(result, decreasing=TRUE);
labels <- yout_test[sidx]; labels <- as.numeric(labels)
TPR = cumsum(labels)/sum(labels); FPR = cumsum(!labels)/sum(!labels);
aTPR[ iter, 1:length(TPR) ] = TPR ; 
aFPR[iter, 1:length(FPR)] = FPR ;
aresults [iter, 1:length(result)] = result[sidx] ;
allpred[,iter] = finalpred;
sen_spe_i <- calculate_results(result, yout_test)
if (length(sen_spe_i) == 0) {stop("Please try a different model (recommended: 'random.forest' or 'svm') ")}
sen_spe_mat_i[iter,] <-  sen_spe_i

}


####################################################################################################################
########################################################## count and give the number of iternations ######################################################################################


cat(paste0(round(iter / n_bootstrap * 100), '% training completed'))
if (iter == n_bootstrap) cat(': Done') else cat('\014')
}


aTPR_pre <- aTPR[, colSums(is.na(aTPR)) != nrow(aTPR)] ## remove NAs
aTPR_final <- colSums(aTPR_pre) / ncol(aTPR_pre) ## taking average

aFPR_pre <- aFPR[, colSums(is.na(aFPR)) != nrow(aFPR)] ## remove NAs
aFPR_final <- colSums(aFPR_pre) / ncol(aFPR_pre) ## taking average

##### plot tpr and fpr to get ROC curve ##### 
df = as.data.frame(cbind(aTPR_final, aFPR_final))
ggplot2::ggplot(df, ggplot2::aes(x=aFPR_final, y=aTPR_final, color="red")) + ggplot2::geom_line()

if (ml_model == 'random.forest') { #taking average of feature important
impvar_mat_avg <- rowSums(impvar_mat)/ncol(impvar_mat) 
top20 <- impvar_mat_avg[order(impvar_mat_avg, decreasing = T)][1:20]

pcorr = matrix( 0, 20,1) ;
for (i in 1:20)
{
pcorr[i] = cor(yout_training ,  ascores_training[ ,names(top20[i])]) ;
}

top20_ti_id_pre <- as.matrix(names(top20))
display_top_pred <- matrix("", 20, 4); 
colnames(display_top_pred) = c('id', 'type', 'name', 'tissue')
display_top_pred[1:20, 1] <- top20_ti_id_pre

for (top_i in 1:nrow(display_top_pred)) {
    replace1 <- which(meta_peak_name_tissue[,1] %in% display_top_pred[top_i,1] )
    if (length(replace1) != 0) {
    display_top_pred[top_i, c(3,4)] <- meta_peak_name_tissue[replace1, c(2,3)]}
    repalce_2 <- which(meta_peak[,1] %in% display_top_pred[top_i,1])
    display_top_pred[top_i, 2] <- meta_peak[repalce_2,2]
}
}

allpred_final <-  rowSums(allpred)/nrow(allpred)
unionPeaks_scores <- cbind(unionPeaks[,4], as.numeric(allpred_final))

i1 <- seq(min(allpred_final), max(allpred_final), by = max(allpred_final) / 10) ##creating intervels for looping

tpr_mat <- matrix(0, length(i1), 2)
for (cut_off in 1:length(i1)) {
     
    num_pred_genes <- length(which(genes %in% unionPeaks[which(allpred_final >= i1[cut_off]),4])) # taking the number of positive genes at the cut off
    tpr = num_pred_genes/ length(which(allpred_final >= i1[cut_off]))
    tpr_mat[cut_off, 1]  = i1[cut_off]
    tpr_mat[cut_off, 2] = tpr

}

tpr_final = tpr_mat[which(tpr_mat[,2] == max(tpr_mat[,2])),1]
predicted_genes_pre_unorder <- unionPeaks_scores[which(allpred_final >= tpr_final),]
predicted_genes_pre_order <- predicted_genes_pre_unorder[order(predicted_genes_pre_unorder[,2], decreasing = T),]
predicted_genes_pre <- unique(predicted_genes_pre_order[,1])
if (length(which(predicted_genes_pre %in% genes)) != 0) {
final_predicted_genes <- predicted_genes_pre[-which(predicted_genes_pre %in% genes)]
} else {final_predicted_genes <- predicted_genes_pre}


sen_spe_mat_final <- colSums(sen_spe_mat_i)/nrow(sen_spe_mat_i) #display and return
print(sen_spe_mat_final)

if (ml_model == 'linear.regression' || ml_model == 'logistic.regression')  {
if(length(final_predicted_genes) > 200) {
    message("Try a different ML model for better result\nrecommended: ml_model = 'xg.boost' or 'svm' ")}
if(length(final_predicted_genes) < 10) {
    message("The number of predicted genes is less");
    message("To fine-tune the predictions:\n Increase the number of input genes\n Or decrease/increase 'n_bootstrap' ");
    message("Try a different ML model\nRecommended: 'random.forest' or 'svm' ") }
if(length(final_predicted_genes) < 150) {message('Predicted genes:'); print(final_predicted_genes)}

}

if (ml_model == 'svm') {
    if(length(final_predicted_genes) > 200) {
        message("Number of predicted genes is too high\nfine-tune the model by increasing/decreasing the n_bootsrap");
        message("Or try a different ML model for better result:\n ml_model = 'xg.boost' or 'random.forest' ")}
    if(length(final_predicted_genes) < 10) {
        message("The number of predicted genes is less");
        message("To fine-tune the predictions:\n Increase the number of input genes\n Or decrease/increase 'n_bootstrap' ");
        message("Try a different ML model: ml_model = 'random.forest' or 'xg.boost' ") }
    if(length(final_predicted_genes) < 150) {message('Predicted genes:'); print(final_predicted_genes)}
}

if (ml_model == 'xg.boost') {
    if(length(final_predicted_genes) > 200) {
        message("Number of predicted genes is too high\nfine-tune the model by increasing/decreasing the n_bootsrap");
        message("Or try a different ML model for better results:\n ml_model =  'svm' or 'random.forest'")}
    if(length(final_predicted_genes) < 10) {
        message("The number of predicted genes is less");
        message("To fine-tune the predictions:\n Increase the number of input genes\n Or decrease/increase 'n_bootstrap' ");
        message("Try a different ML model: ml_model = 'random.forest' or 'svm' ") }
    if(length(final_predicted_genes) < 150) {message('Predicted genes:'); print(final_predicted_genes)}
}

################## return the variables (results) ########## 
if (ml_model == 'random.forest') {
    result_list_rf  = list(final_predicted_genes, sen_spe_mat_final, display_top_pred)
    message('Top predictors:')
    print(display_top_pred)
    if(length(final_predicted_genes) > 200) {
        message("Number of predicted genes is too high\nfine-tune the model by increasing/decreasing the n_bootsrap");
        message("Or try a different ML model for better results:\n ml_model =  'svm' or 'xg.boost'")}
    if(length(final_predicted_genes) < 10) {
        message("The number of predicted genes is less");
        message("To fine-tune the predictions:\n Increase the number of input genes\n Or decrease/increase 'n_bootstrap' ");
        message("Try a different ML model: ml_model = 'svm' or 'xg.boost' ") }
    if(length(final_predicted_genes) < 150 & length(final_predicted_genes) >! 1) {
        message('Predicted genes:');  print(final_predicted_genes) }

    return(result_list_rf) }
        else {return(list(final_predicted_genes, sen_spe_mat_final))}


}
