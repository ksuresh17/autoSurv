# Function that computes survival predictions for a new data set ----------
#' Title
#'
#' @param object
#' @param newdata
#' @param times
#' @param timeVar
#' @param statusVar
#'
#' @return
#' @export
#'
#' @examples
Score_newdata <- function(object, newdata, times, timeVar, statusVar) {
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

    list("auc" = scoreAUC, "brier" = scoreBrier, "pred_probabilities" = preds, "models" = object$models)
}
