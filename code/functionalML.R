## this is the main script, where different machine learning models are trained to make predicton
## Note: We have given commands to run only random forest model, that is, ml_model = 3

args=(commandArgs(TRUE))
library("ggplot2");
library("randomForest"); 
library("glmnet") ; 
library("ROCR")
library("xgboost") ;
library("mltest")
library("SHAPforxgboost")

{

args=(commandArgs(TRUE))
   for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }

}

peakscores <<- matrix(NA, 89747, 0) # reading feature file
for (p_file in 1:42) {
file_name = paste('peakall', p_file, '.rds', sep = "") 
peakscore_path <- paste("../inst/extdata", file_name, sep = "/")
peakscore_n <- readRDS(peakscore_path)
peakscores <<- cbind(peakscores, peakscore_n)
}
pathwaygenes <- read.table("../data/pathway-genes.gmt.csv", sep = "\t", header = F) #reading gene-set file
unionPeaks <- read.table("../data/unionpeak2.0", header=F) #reading file containing curated promoters of genes 
up = as.matrix(unionPeaks) ;
pdim = dim(pathwaygenes) ;
updim = dim(unionPeaks) ;
rsquared <- matrix(0, pdim[1], 1);
allcor = rsquared ;

numfunc = endi - starti +1; #initialing for for loop and creating matricies where results will be saved
allpred = matrix(0, updim[1], numfunc); colnames(allpred) <- pathwaygenes[c(starti:endi), 1]
aTPR = matrix(NA , numfunc , 1000) ;
aFPR = matrix(NA , numfunc , 1000) ;
aresults = matrix(NA , numfunc , 1000) ;
alltfs <- matrix(NA, 20, numfunc) ; 
sen_spe_cut <- matrix(NA, 1, 10) ; colnames(sen_spe_cut) <- c("function-ID", "function_name", "cutoff", "accuracy", "balanced-accuracy", "recall/sensitivity", "specificity", "F1", "MCC", "error-rate") ; 

for (i in starti:endi) # interating through all the functions
 {  tryCatch({
idx = i - starti + 1 ;	 
  genes = pathwaygenes[i, 3:pdim[2]] ; #selecting positive genes
  pos = which(genes != "") ;
  genes = as.matrix(genes[pos]) ; #selecting positive genes
  pos = which( up[,4] %in% genes ) ;
  posScores = peakscores[pos,] ; #selecting the peakscore values (features) for the genes 
  posdim = dim(posScores) ;
#  set.seed(2014)
  negRows = sample(1:updim[1], length(pos), 1) ; #selecting 
  negScores = peakscores[negRows,] ; #selecting the peakscore values (features) for the negative genes 
  ascores = rbind(posScores , negScores) ;
  ascores = as.matrix(ascores)
  yout = matrix(0, (posdim[1] * 2), 1) ;
  yout[1:posdim[1]] =1 ;


samp_yout<- floor(0.75 * nrow(yout))  #dividing the data into training and test
train_yout<-sample(seq_len(nrow(yout)), size = samp_yout); train_yout <- as.matrix(train_yout)
test_yout =  sample((1:nrow(yout))[-train_yout]);
yout_test <- yout[test_yout]; yout_test <- as.numeric(yout_test)
ascores_test <- ascores[test_yout,] ;
yout_training <- yout[train_yout]; yout_training <- as.numeric(yout_training)
ascores_training <- ascores[train_yout,]

#logistic regression model
 if(MLmodel == 1) { #different types of models can be trained depending on the number mentioned in the script
    present0 <- grep('lnr_result', ls(envir=.GlobalEnv))
    if (length(present0) == 0 ) {
    system("mkdir lnr_result")}

 yout_test <- data.matrix(yout_test); yout_test <- as.numeric(yout_test)
 yout_training <- data.matrix(yout_training)
 ascores_training <- data.matrix (ascores_training)
 lm_train1 <- cv.glmnet(ascores_training , yout_training, alpha = 0, family="binomial")
 glm_train <- glmnet(ascores_training, yout_training, alpha = 0, lambda = lm_train1$lambda.min,family = "binomial" )

 newx = data.matrix(ascores_test)
 result = predict( newx= ascores_test, glm_train, s='lambda.min') ;
 finalpred = predict( glm_train, data.matrix(peakscores)) ; #predicting the probablilities for all the promoters of the genes
 
 sidx = order(result, decreasing=TRUE) ;
 labels <- yout_test[sidx] ;

 TPR=cumsum(labels)/sum(labels); FPR= cumsum(!labels)/sum(!labels); #True positive rate and false positive rate calculation
 aTPR[ idx,1:length(TPR) ] = TPR ;
 aFPR[ idx,1:length(FPR) ] = FPR ;
 aresults[ idx,1:length(result) ] = result[sidx] ;
 allpred[,idx] = finalpred ;


 predicate <- prediction(result, yout_test)
roc.perf = performance(predicate, measure = "tpr", x.measure = "fpr")
perf <- performance(predicate,"tpr","fpr")
opt.cut = function(perf, predicate){ #function to calculate best cutoff for highest balanced accuracy
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, predicate@cutoffs)
}
cutoff <- c(idx,opt.cut(roc.perf, predicate ))
predicted_result <- ifelse(result>cutoff[4], 1, 0)
x_report <- ml_test(predicted_result, yout_test) ; funcID <- starti + idx -1
report <- c(funcID, pathwaygenes[i,1], cutoff[4], x_report$accuracy, x_report$balanced.accuracy[2], x_report$recall[2], x_report$specificity[2], x_report$F1[2], x_report$MCC[2], x_report$error.rate)
report <- as.matrix(report) ; report <- t(report)
sen_spe_cut <- rbind(sen_spe_cut, report)
cat(paste0(round(i / numfunc * 100), ' gene-sets completed'))
if (i == numfunc) cat(': Done') else cat('\014')

}


if(MLmodel == 2) { #linerar regression
	present0 <- grep('lgr_result', ls(envir=.GlobalEnv))
    if (length(present0) == 0 ) {
    system("mkdir lgr_result")}

 yout_test <- data.matrix(yout_test)
 yout_training <- data.matrix(yout_training)
 lm_train2 <- cv.glmnet(ascores_training , yout_training, alpha = 1)
 glm_train <- glmnet(ascores_training, yout_training, alpha = 1, lambda = lm_train2$lambda.min)

 newx = data.matrix(ascores_test)
 result = predict( newx= ascores_test, glm_train, s='lambda.min') ;

finalpred = predict( glm_train , data.matrix(peakscores)) ; #predicting the probablilities for all the promoters of the genes
sidx = order(result, decreasing=TRUE) ;
labels <- yout_test[sidx]
TPR=cumsum(labels)/sum(labels); #True positive rate and false positive rate calculation
FPR= cumsum(!labels)/sum(!labels) ;
aTPR[ idx,1:length(TPR) ] = TPR ;
aFPR[ idx,1:length(FPR) ] = FPR ;
aresults[ idx,1:length(result) ] = result[sidx] ;
allpred[,idx] = finalpred ;


predicate <- prediction(result, yout_test)
roc.perf = performance(predicate, measure = "tpr", x.measure = "fpr")
perf <- performance(predicate,"tpr","fpr")
opt.cut = function(perf, predicate){ #function to calculate best cutoff for highest balanced accuracy
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, predicate@cutoffs)
}
cutoff <- c(idx,opt.cut(roc.perf, predicate ))
predicted_result <- ifelse(result>cutoff[4], 1, 0)
x_report <- ml_test(predicted_result, yout_test) ; funcID <- starti + idx -1
report <- c(funcID, pathwaygenes[i,1], cutoff[4], x_report$accuracy, x_report$balanced.accuracy[2], x_report$recall[2], x_report$specificity[2], x_report$F1[2], x_report$MCC[2], x_report$error.rate)
report <- as.matrix(report) ; report <- t(report)
sen_spe_cut <- rbind(sen_spe_cut, report)
cat(paste0(round(i / numfunc * 100), ' gene-sets completed'))
if (i == numfunc) cat(': Done') else cat('\014')

}

if(MLmodel == 3) {
	present0 <- grep('rf_result', system("ls", intern = T))
    if (length(present0) == 0) {
    system("mkdir rf_result")}
    
    set.seed(981998)
lm_train3 <- randomForest(yout_training ~., ntree = 50 , data=data.frame(ascores_training)) ;
predtrain3 = predict( lm_train3, data.frame(ascores_training)) ;

result = predict( lm_train3 , data.frame(ascores_test)) ;
finalpred = predict( lm_train3 , data.frame(peakscores)) ; #predicting the probablilities for all the promoters of the genes
sidx = order(result, decreasing=TRUE) ;
labels <- yout_test[sidx]; labels <- as.numeric(labels)
TPR=cumsum(labels)/sum(labels); FPR= cumsum(!labels)/sum(!labels) ; #True positive rate and false positive rate calculation
aresults[ idx,1:length(result) ] = result[sidx] ; #appending the result of each function to the matrix
aTPR[ idx,1:length(TPR) ] = TPR ; #appending the result of each function to the matrix
aFPR[ idx,1:length(FPR) ] = FPR ; #appending the result of each function to the matrix
allpred[,idx ] = finalpred ; #appending the result of each function to the matrix

imp <- importance(lm_train3, scale = TRUE )  #feature importance 
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
top20 <-impvar[1:20]
alltfs[1:length(top20), idx] = top20 ;  #appending top predictors

predicate <- prediction(result, yout_test)
roc.perf = performance(predicate, measure = "tpr", x.measure = "fpr")
perf <- performance(predicate,"tpr","fpr")
opt.cut = function(perf, predicate){ #function to calculate best cutoff for highest balanced accuracy
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, predicate@cutoffs)
}
cutoff <- c(idx,opt.cut(roc.perf, predicate ))
predicted_result <- ifelse(result>cutoff[4], 1, 0)
x_report <- ml_test(predicted_result, yout_test) ; funcID <- starti + idx -1  ;
report <- c(funcID, pathwaygenes[i,1], cutoff[4], x_report$accuracy, x_report$balanced.accuracy[2], x_report$recall[2], x_report$specificity[2], x_report$F1[2], x_report$MCC[2], x_report$error.rate)
report <- as.matrix(report) ; report <- t(report)
sen_spe_cut <- rbind(sen_spe_cut, report) #gives complete report for the model

cat(paste0(round(i / numfunc * 100), ' gene-sets completed'))
if (i == numfunc) cat(': Done') else cat('\014')


}

 if(MLmodel == 4) {
	present0 <- grep('svm_result', ls(envir=.GlobalEnv))
    if (length(present0) == 0 || length(present1) == 0) {
    system("mkdir svm_result")}

 lm_train4 <- svm(yout_training ~ . , data= data.frame(ascores_training) )
predtrain4 = predict( lm_train4 , data.frame(ascores_training)) ;

result = predict( lm_train4 , data.frame(ascores_test)) ;
finalpred = predict( lm_train4, data.frame(peakscores)) ; #predicting the probablilities for all the promoters of the genes
sidx = order(result, decreasing=TRUE) ;
labels <- yout_test[sidx]
TPR=cumsum(labels)/sum(labels); FPR= cumsum(!labels)/sum(!labels) ; #True positive rate and false positive rate calculation
aTPR[ idx,1:length(TPR) ] = TPR ;
aFPR[ idx,1:length(FPR) ] = FPR ;
aresults[ idx,1:length(result) ] = result[sidx] ;
allpred[,idx ] = finalpred ;

predicate <- prediction(result, yout_test)
roc.perf = performance(predicate, measure = "tpr", x.measure = "fpr")
perf <- performance(predicate,"tpr","fpr")
opt.cut = function(perf, predicate){ #function to calculate best cutoff for highest balanced accuracy
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, predicate@cutoffs)
}
cutoff <- c(idx,opt.cut(roc.perf, predicate ))
predicted_result <- ifelse(result>cutoff[4], 1, 0)
x_report <- ml_test(predicted_result, yout_test) ; funcID <- starti + idx -1 ; 
report <- c(funcID, pathwaygenes[i,1], cutoff[4], x_report$accuracy, x_report$balanced.accuracy[2], x_report$recall[2], x_report$specificity[2], x_report$F1[2], x_report$MCC[2], x_report$error.rate)
report <- as.matrix(report) ; report <- t(report)
sen_spe_cut <- rbind(sen_spe_cut, report)

cat(paste0(round(i / numfunc * 100), ' gene-sets completed'))
if (i == numfunc) cat(': Done') else cat('\014')


 }

if (MLmodel ==5) {
	present0 <- grep('xgb_result', ls(envir=.GlobalEnv))
    if (length(present0) == 0 || length(present1) == 0) {
    system("mkdir xgb_result")}

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

dtrain <- xgb.DMatrix(ascores_training, label = (yout_training))
dvalid <- xgb.DMatrix(data = ascores_validation, label = yout_validation)
dtest <- xgb.DMatrix(ascores_test, label = yout_test)

lowest_error_list = list()
parameters_list = list()


set.seed(20)
for (iter in 1:10){
  param <- list(booster = "gbtree",
                objective = "binary:logistic",
                max_depth = sample(2:8, 1),
                eta = runif(1, .01, .3),
                subsample = runif(1, .7, 1),
                colsample_bytree = runif(1, .8, 1),
                min_child_weight = sample(0:10, 1)
  )
  parameters <- as.data.frame(param)
  parameters_list[[iter]] <- parameters
}

parameters_df = do.call(rbind, parameters_list)

for (row in 1:nrow(parameters_df)){
  set.seed(20)
  mdcv <- xgb.train(data=dtrain,
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
set.seed(20)
params <- list(booster = "gbtree", 
               objective = "binary:logistic",
               max_depth = randomsearch[1,]$max_depth,
               eta = randomsearch[1,]$eta,
               subsample = randomsearch[1,]$subsample,
               colsample_bytree = randomsearch[1,]$colsample_bytree,
               min_child_weight = randomsearch[1,]$min_child_weight)
lm_train5 <- xgb.train(params = params,
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
finalpred <- predict(lm_train5, as.matrix(peakscores)); #predicting the probablilities for all the promoters of the genes
sidx <- order(result, decreasing=TRUE);
labels <- yout_test[sidx]; labels <- as.numeric(labels)
TPR = cumsum(labels)/sum(labels); FPR = cumsum(!labels)/sum(!labels); #True positive rate and false positive rate calculation
aTPR[ idx, 1:length(TPR) ] = TPR ; 
aFPR[idx, 1:length(FPR)] = FPR ;
aresults [idx, 1:length(result)] = result[sidx] ;
allpred[,idx] = finalpred;

imp <- xgb.importance(feature_names=NULL, model = lm_train5)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
top20 <-impvar[1:20]
alltfs[1:length(top20), idx] = top20 ;

predicate <- prediction(result, yout_test)
roc.perf = performance(predicate, measure = "tpr", x.measure = "fpr")
perf <- performance(predicate,"tpr","fpr")
opt.cut = function(perf, predicate){ #function to calculate best cutoff for highest balanced accuracy
    cut.ind = mapply(FUN=function(x, y, p){
        d = (x - 0)^2 + (y-1)^2
        ind = which(d == min(d))
        c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
            cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, predicate@cutoffs)
}
cutoff <- c(idx,opt.cut(roc.perf, predicate ))
predicted_result <- ifelse(result>cutoff[4], 1, 0)
x_report <- ml_test(predicted_result, yout_test) ; funcID <- starti + idx -1 ;
report <- c(funcID, pathwaygenes[i,1], cutoff[4], x_report$accuracy, x_report$balanced.accuracy[2], x_report$recall[2], x_report$specificity[2], x_report$F1[2], x_report$MCC[2], x_report$error.rate)
report <- as.matrix(report) ; report <- t(report)
sen_spe_cut <- rbind(sen_spe_cut, report)

shap_values <- shap.values(xgb_model = lm_train5, X_train = ascores_training)
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = ascores_training)
variable <- unique(shap_long$variable)
top10_shap <- t(as.matrix(variable[1:10]))
top10_shap <- cbind(i,top10_shap)

shap_matrix <- rbind(shap_matrix, top10_shap)



cat(paste0(round(i / numfunc * 100), ' gene-sets completed'))
if (i == numfunc) cat(': Done') else cat('\014')

}


allcor[i] = cor(result, yout_test) ;
},
 error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}
sen_spe_cut <- sen_spe_cut[-1,]
shap_matrix <- shap_matrix[-1,]

write.table( aTPR, file=TPRf) ;
write.table(  aFPR , file=FDRf) ;
write.table(allcor, file= outf ) ;
write.table(allpred, file=predf ) ;
if (MLmodel == 5 | MLmodel == 3) {
write.table(alltfs, file=toptf);
}
write.table(sen_spe_cut, file=cfmtx);
