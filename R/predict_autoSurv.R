#' Discrete-time survival prediction models 
#' 
#' @description 
#' Compute the predicted survival probabilities for trained discrete-time models
#'
#' @section Authors:
#' Krithika Suresh (\email{krihtika.suresh@@cuanschutz.edu})
#'
#' @param object autoSurv() object that contains trained models 
#' @param newdata data set on which you want to obtain survival predictions 
#' @param times vector of time points at which survival probability is predicted
#' @param timeVar string corresponding to the variable name of the time-to-event outcome. Used to compute performance metrics
#' @param statusVar string corresponding to the variable name of the status indicator. Used to compute performance metrics
#'
#' @return A list with predicted probabilities and performance metrics
#' @export

predict_autoSurv <- function(object, newdata, times, timeVar, statusVar) {
    # create survival formula for censoring
    .form_cens <- as.formula(paste("Surv(", timeVar,",", statusVar, ") ~ 1"))
    #list of models that were fit to the original data
    trainModels <- names(object$models)
    #store predictions for assessment
    preds <- list()
    
    if("cox" %in% trainModels) {
        preds_cox <-  matrix(NA, nrow(newdata), length(times))
        for (i in 1:length(times)){
            preds_cox[,i] <- predictCox(object$models$cox, dat=newdata, w=times[i])
        }
        preds[["cox"]] <- preds_cox
    }
    
    if("rsf" %in% trainModels) {
        preds_rsf <-  matrix(NA, nrow(newdata), length(times))
        for (i in 1:length(times)){
            preds_rsf[,i] <- predictRSF(object$models$rsf, newdata=newdata, times=times[i])
        }
        preds[["rsf"]] <- preds_rsf
    }
    
    if ("glm" %in% trainModels){
        preds_glm <- predBinSurv(testDat = newdata, delta.upper = object$models$glm[["delta.upper"]],
                                 model = object$models$glm[["model"]], times)
        preds[["glm"]] <- preds_glm
    }
    
    if ("gbm" %in% trainModels){
        preds_gbm <- predBinSurv(testDat = newdata, delta.upper = object$models$gbm[["delta.upper"]],
                                 model = object$models$gbm[["model"]], times)
        preds[["gbm"]] <- preds_gbm
    }
    
    if ("gam" %in% trainModels){
        preds_gam <- predBinSurv(testDat = newdata, delta.upper = object$models$gam[["delta.upper"]],
                                 model = object$models$gam[["model"]], times)
        preds[["gam"]] <- preds_gam
    }
    
    if ("glmnet" %in% trainModels){
        preds_glmnet <- predBinSurv(testDat = newdata, delta.upper = object$models$glmnet[["delta.upper"]],
                                    model = object$models$glmnet[["model"]], times)
        preds[["glmnet"]] <- preds_glmnet
    }
    
    if ("svm" %in% trainModels){
        preds_svm <- predBinSurv(testDat = newdata, delta.upper = object$models$svm[["delta.upper"]],
                                 model = object$models$svm[["model"]], times)
        preds[["svm"]] <- preds_svm
    }
    
    if ("cforest" %in% trainModels){
        preds_cforest <- predBinSurv(testDat = newdata, delta.upper = object$models$cforest[["delta.upper"]],
                                     model = object$models$cforest[["model"]], times)
        preds[["cforest"]] <- preds_cforest
    }
    
    if ("nnet" %in% trainModels){
        preds_nnet <- predBinSurv(testDat = newdata, delta.upper = object$models$nnet[["delta.upper"]],
                                  model = object$models$nnet[["model"]], times)
        preds[["nnet"]] <- preds_nnet
    }
    
    temp_Score <- riskRegression::Score(preds,
                                        formula=.form_cens,
                                        data=newdata,
                                        times=times,
                                        summary="ibs")
    
    scoreAUC <- temp_Score$AUC$score
    scoreBrier <- temp_Score$Brier$score
    scoreBrier$R2 <- 1-scoreBrier$Brier/scoreBrier$Brier[which(scoreBrier$model=="Null model")]
    scoreBrier$R2_IBS <- 1-scoreBrier$IBS/scoreBrier$IBS[which(scoreBrier$model=="Null model")]
    
    scoreALL <- data.frame(model = as.character(scoreAUC$model),
                           AUC = scoreAUC$AUC,
                           Brier = scoreBrier$Brier[-which(scoreBrier$model=="Null model")],
                           R2 = scoreBrier$R2[-which(scoreBrier$model=="Null model")],
                           IBS = scoreBrier$IBS[-which(scoreBrier$model=="Null model")],
                           R2_IBS = scoreBrier$R2_IBS[-which(scoreBrier$model=="Null model")])
    
    list("auc" = scoreAUC,
         "brier" = scoreBrier,
         "metrics" = scoreALL,
         "pred_probabilities" = lapply(preds, function(x) 1-x),
         "models" = object$models)
}
