### TODO ###
# Extract models from discrete time. may need to rework how function works
# figure out why intervals > 10 have issues. Krithika looking at this.
# fix id, fold names to be unique
# fix results table to be presentable
# fix warning messages during training
# add neural net and SVM
# diagnose issues with glmnet poor performance.
# add hyperparameter tuning for other discrete time model, gbm is done
# maybe add cross validation
# clean code for presentation on github

#' autoSurv
#'
#' @param timeVar string corresponding to the variable name of the time-to-event outcome
#' @param statusVar string corresponding to the variable name of the status indicator
#' @param data data frame containing covariates for prediction
#' @param times vector of time points at which probably of event is predicted
#' @param trainModels list of models to train. options: "cox","rsf","glm","gam","gbm","glmnet","svm","cforest","nnet"
#' @param bins number of time bins
#' @param cens string indicating which method to use for handling censored observations. options: "same","prev","half"
#' @param testProp proportion of observations to be randomly selected for the test set
#' @param seed integer random seed for reproducing results
#' @param init number of randomly chosen points to sample the target function before Bayesian Optimization fitting the Gaussian Process
#' @param iter total number of times the Bayesian Optimization is to be repeated
#' @param cv enable or disable cross-validation
#' @param cvFold number of folds to use in k-fold cross validation
#' @param verbose.opt enable or disable printing of the Bayesian optimization progress
#'
#' @return
#' @export
#'
#' @examples
autoSurv <- function(timeVar,
                     statusVar,
                     data,
                     times = NULL,
                     trainModels = c("cox","rsf","glm","gbm","glmnet","svm","cforest","nnet"),
                     bins.lower = 5,
                     bins.upper = 10,
                     cens = "same",
                     testProp = 0.2,
                     seed = 1,
                     init = 10,
                     iter = 20,
                     cv = TRUE,
                     cvFold = 5,
                     verbose.opt = TRUE
) {

    # split into train / test
    if (testProp < 1 & testProp > 0){
        set.seed(seed)
        testIndexes <- sample(x = 1:nrow(data),size = nrow(data)*testProp)
        testDat <- data[testIndexes, ]
        trainDat <- data[-testIndexes, ]
    } else if (testProp == 0){
        warning("testProp is set to 0. Performance results will be measured
                using the same data as training. Note: These results may be
                overfitted")
        testDat <- data
        trainDat <- data
    } else {
        stop("testProp must be a non-negative value less than 1")
    }

    # Assign folds for internal cross validation
    if(cv){
        if(cvFold > 1){
            set.seed(seed)
            internal_fold <- sample(seq_len(cvFold),nrow(trainDat),1/cvFold)
        } else{
            stop("cvFold should be a positive integer greater than 1.")
        }
    } else {
        cvFold <- 1
        if(testProp > 0 & testProp < 1){
            set.seed(seed)
            testIndexes <- sample(x = 1:nrow(trainDat),size = nrow(trainDat)*testProp)
            fold_test <- trainDat[testIndexes,]
            fold_train <- trainDat[-testIndexes,]
        } else {
            fold_test <- data
            fold_train <- data
        }

    }

    #if prediction time points are missing use the median of event times
    if(is.null(times)){
        times <- quantile(trainDat[[timeVar]][trainDat[[statusVar]]==1], 0.5, names=FALSE)
    }

    #if more than one time specified, then use last time for optimizing
    if(length(times)>1) {
        w_pred <- tail(times,1)
    } else {
        w_pred <- times
    }

    # create survival formula for continuous survival models
    .form <- as.formula(paste("Surv(", timeVar, ",", statusVar, ") ~ ."))
    # create survival formula for censoring
    .form_cens <- as.formula(paste("Surv(", timeVar,",", statusVar, ") ~ 1"))

    #initialize list for saving models and predictions
    saveModels <- list()
    optParams <- list()
    preds <- list()

    # cox model
    if ("cox" %in% trainModels){
        coxphModel <- survival::coxph(.form,
                                      data=trainDat)

        pred_cox <-  matrix(NA, nrow(testDat), length(times))
        for (i in 1:length(times)){
            pred_cox[,i] <- predictCox(coxphModel, dat=testDat, w=times[i])
        }
        preds[["cox"]] <- pred_cox
        saveModels$cox <- coxphModel
    }

    # Random Survival Forests
    if ("rsf" %in% trainModels){

        rf_cv <- function(nodesize, mtry){
            fold_briers <- NULL
            for (i in seq_len(cvFold)){
                if(cv){
                    fold_test <- trainDat[internal_fold == i,]
                    fold_train <- trainDat[internal_fold != i,]
                }


                rfsrcModel <- randomForestSRC::rfsrc(.form,
                                                     nodesize = nodesize,
                                                     mtry = mtry,
                                                     data=fold_train)
                brier <- NULL
                temp_Score <- riskRegression::Score(list( "RSF" = rfsrcModel),
                                                    formula=.form_cens,
                                                    data=fold_test,
                                                    times=times,
                                                    summary="ibs")
                if(length(times)>1) {
                    brier <- tail(temp_Score$Brier$score$IBS[temp_Score$Brier$score$model == "RSF"], 1)
                } else {
                    brier <- temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "RSF"]
                }

                fold_briers <- c(fold_briers, brier)
            }

            neg_brier <- -mean(fold_briers)

            list(Score = neg_brier,
                 Pred = 0)
        }

        print("Optimizing RSF")

        OPT_Res <- rBayesianOptimization::BayesianOptimization(rf_cv,
                                                               bounds = list(nodesize = c(1L, 15L),
                                                                             mtry = c(2L, as.integer(ncol(trainDat)-2))),
                                                               init_grid_dt = NULL,
                                                               init_points = init,
                                                               n_iter = iter,
                                                               acq = "ucb",
                                                               kappa = 2.576,
                                                               eps = 0.0,
                                                               verbose = verbose.opt)

        rfsrcModel <- randomForestSRC::rfsrc(.form,
                                             nodesize = OPT_Res$Best_Par[["nodesize"]],
                                             mtry = OPT_Res$Best_Par[["mtry"]],
                                             data=trainDat)

        pred_rsf <-  matrix(NA, nrow(testDat), length(times))
        for (i in 1:length(times)){
            pred_rsf[,i] <- predictRSF(rfsrcModel, newdata=testDat, times=times[i])
        }
        preds[["rsf"]] <- pred_rsf
        saveModels$rsf <- rfsrcModel
        optParams$rsf <- OPT_Res$Best_Par
    }

    #Discrete time models

    # GLM discrete time model
    if ("glm" %in% trainModels){
        glm_cv <- function(intervals){
            fold_briers <- NULL
            for (i in seq_len(cvFold)){
                if(cv){
                    fold_test <- trainDat[internal_fold == i,]
                    fold_train <- trainDat[internal_fold != i,]
                }

                pred.GLM <- glmBinModel(trainDat = fold_train,
                                        timeVar = timeVar,
                                        statusVar = statusVar,
                                        w=w_pred,
                                        bins=intervals,
                                        cens=cens)

                pred.GLM.mat <- predBinSurv(testDat = fold_test,
                                            delta.upper = pred.GLM[["delta.upper"]],
                                            model = pred.GLM[["model"]],
                                            times)

                temp_Score <- riskRegression::Score(list("GLM" = pred.GLM.mat),
                                                    formula = .form_cens,
                                                    data = fold_test,
                                                    times = times,
                                                    summary = "ibs")

                if(length(times)>1) {
                    brier <- tail(temp_Score$Brier$score$IBS[temp_Score$Brier$score$model == "GLM"], 1)
                } else {
                    brier <- temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "GLM"]
                }
                fold_briers <- c(fold_briers, brier)
            }

            neg_brier <- -mean(fold_briers)

            list(Score = neg_brier,
                 Pred = 0)
        }

        print("Optimizing DiscreteTime-GLM")

        OPT_Res_GLM <- rBayesianOptimization::BayesianOptimization(glm_cv,
                                                                   bounds = list(intervals = c(as.integer(bins.lower), as.integer(bins.upper))),
                                                                   init_grid_dt = NULL,
                                                                   init_points = init,
                                                                   n_iter = iter,
                                                                   acq = "ucb",
                                                                   kappa = 2.576,
                                                                   eps = 0.0,
                                                                   verbose = verbose.opt)

        OPT_pred.GLM <- glmBinModel(trainDat = trainDat,
                                    timeVar = timeVar,
                                    statusVar = statusVar,
                                    w=w_pred,
                                    bins=OPT_Res_GLM$Best_Par[["intervals"]],
                                    cens=cens)

        preds_GLM <- predBinSurv(testDat = testDat,
                                 delta.upper = OPT_pred.GLM[["delta.upper"]],
                                 model = OPT_pred.GLM[["model"]],
                                 times)

        preds[["glm"]] <- preds_GLM
        saveModels$glm <- OPT_pred.GLM
        optParams$glm <- OPT_Res_GLM$Best_Par
    }

    # GBM discrete time model
    if ("gbm" %in% trainModels){
        gbm_cv <- function(intervals, n.trees, interaction.depth, shrinkage, n.minobsinnode){
            fold_briers <- NULL
            for (i in seq_len(cvFold)){
                if(cv){
                    fold_test <- trainDat[internal_fold == i,]
                    fold_train <- trainDat[internal_fold != i,]
                }

                pred.GBM <- gbmBinModel(trainDat = fold_train,
                                        timeVar = timeVar,
                                        statusVar= statusVar,
                                        w = w_pred,
                                        bins = intervals,
                                        cens = cens,
                                        n.trees = n.trees,
                                        interaction.depth = interaction.depth,
                                        shrinkage = shrinkage,
                                        n.minobsinnode = n.minobsinnode
                )

                pred.GBM.mat <- predBinSurv(testDat = fold_test,
                                            delta.upper = pred.GBM[["delta.upper"]],
                                            model = pred.GBM[["model"]],
                                            times)

                temp_Score <- riskRegression::Score(list("GBM" = pred.GBM.mat),
                                                    formula=.form_cens,
                                                    data=fold_test,
                                                    times=times,
                                                    summary="ibs")

                if(length(times)>1) {
                    brier <- tail(temp_Score$Brier$score$IBS[temp_Score$Brier$score$model == "GBM"], 1)
                } else {
                    brier <- temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "GBM"]
                }

                fold_briers <- c(fold_briers, brier)
            }

            neg_brier <- -mean(fold_briers)

            list(Score = neg_brier,
                 Pred = 0)
        }

        print("Optimizing DiscreteTime-GBM")

        OPT_Res_GBM <- rBayesianOptimization::BayesianOptimization(gbm_cv,
                                                                   bounds = list(intervals = c(as.integer(bins.lower), as.integer(bins.upper)),
                                                                                 n.trees = c(50L, 500L),
                                                                                 interaction.depth = c(1L,3L),
                                                                                 shrinkage = c(0.001,0.1),
                                                                                 n.minobsinnode = c(1L,20L)
                                                                   ),
                                                                   init_grid_dt = NULL,
                                                                   init_points = init,
                                                                   n_iter = iter,
                                                                   acq = "ucb",
                                                                   kappa = 2.576,
                                                                   eps = 0.0,
                                                                   verbose = verbose.opt)

        OPT_pred.GBM <- gbmBinModel(trainDat = trainDat,
                                    timeVar = timeVar,
                                    statusVar = statusVar,
                                    w = w_pred,
                                    bins=OPT_Res_GBM$Best_Par[["intervals"]],
                                    cens=cens,
                                    n.trees = OPT_Res_GBM$Best_Par[["n.trees"]],
                                    interaction.depth = OPT_Res_GBM$Best_Par[["interaction.depth"]],
                                    shrinkage = OPT_Res_GBM$Best_Par[["shrinkage"]],
                                    n.minobsinnode = OPT_Res_GBM$Best_Par[["n.minobsinnode"]]
        )

        preds_GBM <- predBinSurv(testDat = testDat, delta.upper = OPT_pred.GBM[["delta.upper"]],
                                 model = OPT_pred.GBM[["model"]], times)

        preds[["gbm"]] <- preds_GBM
        saveModels$gbm <- OPT_pred.GBM
        optParams$gbm <- OPT_Res_GBM$Best_Par
    }

    # glmnet discrete time model
    if ("glmnet" %in% trainModels){
        glmnet_cv <- function(intervals, alpha, lambda){
            fold_briers <- NULL
            for (i in seq_len(cvFold)){
                if(cv){
                    fold_test <- trainDat[internal_fold == i,]
                    fold_train <- trainDat[internal_fold != i,]
                }

                pred.glmnet <- glmnetBinModel(trainDat = fold_train,
                                              timeVar = timeVar,
                                              statusVar = statusVar,
                                              w = w_pred,
                                              bins = intervals,
                                              cens = cens,
                                              alpha = alpha,
                                              lambda = lambda
                )

                pred.glmnet.mat <- predBinSurv(testDat = fold_test,
                                               delta.upper = pred.glmnet[["delta.upper"]],
                                               model = pred.glmnet[["model"]], times)

                temp_Score <- riskRegression::Score(list("glmnet" = pred.glmnet.mat),
                                                    formula=.form_cens,
                                                    data=fold_test,
                                                    times=times,
                                                    summary="ibs")

                if(length(times)>1) {
                    brier <- tail(temp_Score$Brier$score$IBS[temp_Score$Brier$score$model == "glmnet"], 1)
                } else {
                    brier <- temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "glmnet"]
                }
                fold_briers <- c(fold_briers, brier)
            }

            neg_brier <- -mean(fold_briers)

            list(Score = neg_brier,
                 Pred = 0)
        }

        print("Optimizing DiscreteTime-glmnet")

        OPT_Res_glmnet<- rBayesianOptimization::BayesianOptimization(glmnet_cv,
                                                                     bounds = list(intervals = c(as.integer(bins.lower), as.integer(bins.upper)),
                                                                                   alpha = c(0, 0.5),
                                                                                   lambda = c(0, 0.5)
                                                                     ),
                                                                     init_grid_dt = NULL,
                                                                     init_points = init,
                                                                     n_iter = iter,
                                                                     acq = "ucb",
                                                                     kappa = 2.576,
                                                                     eps = 0.0,
                                                                     verbose = verbose.opt)

        OPT_pred.glmnet <- glmnetBinModel(trainDat = trainDat,
                                          timeVar = timeVar,
                                          statusVar= statusVar,
                                          w = w_pred,
                                          bins=OPT_Res_glmnet$Best_Par[["intervals"]],
                                          cens = cens,
                                          alpha=OPT_Res_glmnet$Best_Par[["alpha"]],
                                          lambda=OPT_Res_glmnet$Best_Par[["lambda"]]
        )

        preds_glmnet <- predBinSurv(testDat = testDat, delta.upper = OPT_pred.glmnet[["delta.upper"]],
                                    model = OPT_pred.glmnet[["model"]], times)

        preds[["glmnet"]] <- preds_glmnet
        saveModels$glmnet <- OPT_pred.glmnet
        optParams$glmnet <- OPT_Res_glmnet$Best_Par
    }

    # svm discrete time model
    if ("svm" %in% trainModels){
        svm_cv <- function(intervals, C, sigma){
            fold_briers <- NULL
            for (i in seq_len(cvFold)){
                if(cv){
                    fold_test <- trainDat[internal_fold == i,]
                    fold_train <- trainDat[internal_fold != i,]
                }

                pred.svm <- svmBinModel(trainDat = fold_train,
                                        timeVar = timeVar,
                                        statusVar = statusVar,
                                        w = w_pred,
                                        bins = intervals,
                                        cens = cens,
                                        C = C,
                                        sigma = sigma
                )

                pred.svm.mat <- predBinSurv(testDat = fold_test,
                                            delta.upper = pred.svm[["delta.upper"]],
                                            model = pred.svm[["model"]], times)

                temp_Score <- riskRegression::Score(list("svm" = pred.svm.mat),
                                                    formula=.form_cens,
                                                    data=fold_test,
                                                    times=times,
                                                    summary="ibs")

                if(length(times)>1) {
                    brier <- tail(temp_Score$Brier$score$IBS[temp_Score$Brier$score$model == "svm"], 1)
                } else {
                    brier <- temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "svm"]
                }
                fold_briers <- c(fold_briers, brier)
            }

            neg_brier <- -mean(fold_briers)

            list(Score = neg_brier,
                 Pred = 0)
        }

        print("Optimizing DiscreteTime-svm")

        OPT_Res_svm <- rBayesianOptimization::BayesianOptimization(svm_cv,
                                                                   bounds = list(intervals = c(as.integer(bins.lower), as.integer(bins.upper)),
                                                                                 C = c(0, 30),
                                                                                 sigma = c(0.01, 0.2)
                                                                   ),
                                                                   init_grid_dt = NULL,
                                                                   init_points = init,
                                                                   n_iter = iter,
                                                                   acq = "ucb",
                                                                   kappa = 2.576,
                                                                   eps = 0.0,
                                                                   verbose = verbose.opt)

        OPT_pred.svm <- svmBinModel(trainDat = trainDat,
                                    timeVar = timeVar,
                                    statusVar= statusVar,
                                    w = w_pred,
                                    bins = OPT_Res_svm$Best_Par[["intervals"]],
                                    cens = cens,
                                    C= OPT_Res_svm$Best_Par[["C"]],
                                    sigma = OPT_Res_svm$Best_Par[["sigma"]]
        )

        preds_svm <- predBinSurv(testDat = testDat, delta.upper = OPT_pred.svm[["delta.upper"]],
                                 model = OPT_pred.svm[["model"]], times)

        preds[["svm"]] <- preds_svm
        saveModels$svm <- OPT_pred.svm
        optParams$svm <- OPT_Res_svm$Best_Par
    }

    # cforest discrete time model
    if ("cforest" %in% trainModels){
        cforest_cv <- function(intervals, mtry){
            fold_briers <- NULL
            for (i in seq_len(cvFold)){
                if(cv){
                    fold_test <- trainDat[internal_fold == i,]
                    fold_train <- trainDat[internal_fold != i,]
                }

                pred.cforest <- cforestBinModel(trainDat = fold_train,
                                                timeVar = timeVar,
                                                statusVar= statusVar,
                                                w = w_pred,
                                                bins = intervals,
                                                cens = cens,
                                                mtry = mtry
                )

                pred.cforest.mat <- predBinSurv(testDat = fold_test,
                                                delta.upper = pred.cforest[["delta.upper"]],
                                                model = pred.cforest[["model"]], times)

                temp_Score <- riskRegression::Score(list("cforest" = pred.cforest.mat),
                                                    formula=.form_cens,
                                                    data=fold_test,
                                                    times=times,
                                                    summary="ibs")

                if(length(times)>1) {
                    brier <- tail(temp_Score$Brier$score$IBS[temp_Score$Brier$score$model == "cforest"], 1)
                } else {
                    brier <- temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "cforest"]
                }
                fold_briers <- c(fold_briers, brier)
            }

            neg_brier <- -mean(fold_briers)

            list(Score = neg_brier,
                 Pred = 0)
        }

        print("Optimizing DiscreteTime-cforest")

        OPT_Res_cforest <- rBayesianOptimization::BayesianOptimization(cforest_cv,
                                                                       bounds = list(intervals = c(as.integer(bins.lower), as.integer(bins.upper)),
                                                                                     mtry = c(2L, as.integer(ncol(trainDat)-2))
                                                                       ),
                                                                       init_grid_dt = NULL,
                                                                       init_points = init,
                                                                       n_iter = iter,
                                                                       acq = "ucb",
                                                                       kappa = 2.576,
                                                                       eps = 0.0,
                                                                       verbose = verbose.opt)

        OPT_pred.cforest <- cforestBinModel(trainDat = trainDat,
                                            timeVar = timeVar,
                                            statusVar = statusVar,
                                            w = w_pred,
                                            bins = OPT_Res_cforest$Best_Par[["intervals"]],
                                            cens = cens,
                                            mtry = OPT_Res_cforest$Best_Par[["mtry"]]
        )

        preds_cforest <- predBinSurv(testDat = testDat,
                                     delta.upper = OPT_pred.cforest[["delta.upper"]],
                                     model = OPT_pred.cforest[["model"]],
                                     times)

        preds[["cforest"]] <- preds_cforest
        saveModels$cforest <- OPT_pred.cforest
        optParams$cforest <- OPT_Res_cforest$Best_Par
    }

    # nnet discrete time model
    if ("nnet" %in% trainModels){
        nnet_cv <- function(intervals, size, decay){

            fold_briers <- NULL
            for (i in seq_len(cvFold)){
                if(cv){
                    fold_test <- trainDat[internal_fold == i,]
                    fold_train <- trainDat[internal_fold != i,]
                }

                pred.nnet <- nnetBinModel(trainDat = fold_train,
                                          timeVar = timeVar,
                                          statusVar= statusVar,
                                          w = w_pred,
                                          bins = intervals,
                                          cens = cens,
                                          size = size,
                                          decay = decay
                )

                pred.nnet.mat <- predBinSurv(testDat = fold_test,
                                             delta.upper = pred.nnet[["delta.upper"]],
                                             model = pred.nnet[["model"]],
                                             times)

                temp_Score <- riskRegression::Score(list("nnet" = pred.nnet.mat),
                                                    formula=.form_cens,
                                                    data=fold_test,
                                                    times=times,
                                                    summary="ibs")

                if(length(times)>1) {
                    brier <- tail(temp_Score$Brier$score$IBS[temp_Score$Brier$score$model == "nnet"], 1)
                } else {
                    brier <- temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "nnet"]
                }
                fold_briers <- c(fold_briers, brier)
            }

            neg_brier <- -mean(fold_briers)

            list(Score = neg_brier,
                 Pred = 0)
        }

        print("Optimizing DiscreteTime-nnet")

        OPT_Res_nnet <- rBayesianOptimization::BayesianOptimization(nnet_cv,
                                                                    bounds = list(intervals = c(as.integer(bins.lower), as.integer(bins.upper)),
                                                                                  size = c(1L, 10L),
                                                                                  decay = c(0.01, 0.5)
                                                                    ),
                                                                    init_grid_dt = NULL,
                                                                    init_points = init,
                                                                    n_iter = iter,
                                                                    acq = "ucb",
                                                                    kappa = 2.576,
                                                                    eps = 0.0,
                                                                    verbose = verbose.opt)

        OPT_pred.nnet <- nnetBinModel(trainDat = trainDat,
                                      timeVar = timeVar,
                                      statusVar = statusVar,
                                      w = w_pred,
                                      bins=OPT_Res_nnet$Best_Par[["intervals"]],
                                      cens=cens,
                                      size = OPT_Res_nnet$Best_Par[["size"]],
                                      decay = OPT_Res_nnet$Best_Par[["decay"]]
        )

        preds_nnet <- predBinSurv(testDat = testDat, delta.upper = OPT_pred.nnet[["delta.upper"]],
                                  model = OPT_pred.nnet[["model"]], times)

        preds[["nnet"]] <- preds_nnet
        saveModels$nnet <- OPT_pred.nnet
        optParams$nnet <- OPT_Res_nnet$Best_Par
    }

    temp_Score <- riskRegression::Score(preds,
                                        formula=.form_cens,
                                        data=testDat,
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
         "models" = saveModels,
         "tuned_params" = optParams
    )
}

# Function that computes survival predictions for a new data set ----------
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


# Function that identifies the endpoints for the discrete intervals --------
createDiscreteIntervals <- function(time, event, bins) {
    #Number of intervals to create (minimum of specified bins and number of events)
    n.delta <- min(bins, length(unique(time[event==1])))
    #Get percentiles for each of the bins (from 0 to 1)
    probs.delta <- seq(from=0, to=1, length.out=(n.delta+1))
    #Get event times for each of the percentiles (0: first event time in dat set, 1: last event time in data set)
    delta.upper <- quantile(time[event==1], probs=probs.delta, names=FALSE)
    #Drop the first interval
    return(delta.upper[-1])
}


# Function that returns binary process data set ---------------------------
# This function was adapted from code from the Github account of Eric Polley
# https://github.com/ecpolley/SuperLearner_Old/blob/master/R/createDiscrete.R
createDiscreteDat <- function(timeInternal, eventInternal, dataX, delta.upper, cens="same") {
    #Assumes t0=0
    delta.lower <- c(0, delta.upper[-length(delta.upper)])
    n.delta <- length(delta.upper)
    IDInternal <- 1:nrow(dataX)
    dat_i <- cbind(IDInternal, timeInternal, eventInternal, dataX)
    interval <- rep(1:length(delta.upper), times=nrow(dataX))

    #Create a long data set that repeats each persons data n.delta (number of intervals) times
    long.dat <- dat_i[rep(seq_len(nrow(dat_i)), each = n.delta), ]

    N.delta <- rep(NA, nrow(long.dat))
    long.dat <- cbind(long.dat, delta.lower, delta.upper, N.delta, interval)

    # Treatment of censored observations
    if(cens=="same") {
        # Include censored observations in the interval in which they were censored
        long.dat$N.delta <- ifelse(long.dat$timeInternal > long.dat$delta.upper, 0,
                                   ifelse(long.dat$eventInternal==1, ifelse(long.dat$timeInternal <= long.dat$delta.lower, NA, 1),
                                          ifelse(long.dat$timeInternal>long.dat$delta.lower, 0, NA)))
    } else if(cens=="prev") {
        # Include censored observations in the previous interval
        long.dat$N.delta <- ifelse(long.dat$timeInternal > long.dat$delta.upper, 0,
                                   ifelse(long.dat$eventInternal==1, ifelse(long.dat$timeInternal <= long.dat$delta.lower, NA, 1),
                                          NA))
    } else if(cens=="half") {
        # Include censored observations in interval if they have survived at least half of that interval, otherwise include in previous interval
        long.dat$N.delta <- ifelse(long.dat$timeInternal > long.dat$delta.upper, 0,
                                   ifelse(long.dat$eventInternal==1, ifelse(long.dat$timeInternal <= long.dat$delta.lower, NA, 1),
                                          ifelse(long.dat$timeInternal>=0.5*(long.dat$delta.upper+long.dat$delta.lower), 0, NA)))
    }

    m <- delta.upper[n.delta]
    # For event time corresponding to last interval include as an event in that interval
    long.dat$N.delta <- ifelse(long.dat$timeInternal == m & long.dat$delta.upper == m,
                               ifelse(is.na(long.dat$N.delta), 0, long.dat$N.delta), long.dat$N.delta)

    # Drops the intervals at which the individual no longer contributes
    long.dat <- long.dat[!is.na(long.dat$N.delta), ]
    # Sort the data
    long.dat <- long.dat[order(long.dat$IDInternal, long.dat$delta.lower),]
    # Set interval as a factor
    long.dat$interval <- as.factor(long.dat$interval)
    # Set outcome as a factor
    long.dat$N.delta <- as.factor(long.dat$N.delta)
    return(long.dat)
}

# Function to create person-time data set ---------------------------------
# dat: (data.frame) data set to convert (with columns named timeVar, statusVar)
# bins: (integer) number of discrete intervals
# w: (integer) prediction interval
genPTdat <- function(dat, timeVar = timeVar, statusVar = statusVar, bins = 5, w, cens) {
    time.ind <- which(eval(timeVar) == colnames(dat))
    event.ind <- which(eval(statusVar) == colnames(dat))

    dat_i <- dat
    dat_i[dat_i$survTime > w, event.ind] <- 0 #ignore events after the prediction window
    dat_i[dat_i$survTime > w, time.ind] <- w #administratively censor at w

    # Identify the times t1, t2,... tJ (assumes t0=0)
    # Don't need to include prediction horizon as a final interval if there are no event times in that interval (prediction close to 0)
    delta.upper <- unique(createDiscreteIntervals(time = dat_i[,time.ind],
                                                  event = dat_i[event.ind],
                                                  bins = bins))
    delta.upper.new <- c(delta.upper[-length(delta.upper)],w)
    delta.upper <- delta.upper.new

    # Creates a new data set where each person has a row corresponding to the discrete time intervals
    # 0: in intervals where alive but does not have event
    # 1: in intervals where experience the event
    # delta.upper and delta.lower columns: are the same for everyone and are the discretized time intervals based on quantiles
    dat_i.X <- createDiscreteDat(timeInternal = dat_i[,time.ind],
                                 eventInternal = dat_i[,event.ind],
                                 dataX = dat_i[,-c(time.ind,event.ind)],
                                 delta.upper = delta.upper,
                                 cens=cens)

    return(dat_i.X)
}

# Function that formats data we will make predictions on ------------------
# Similar to "genPTdat" function but does not create outcome variable and creates
# a row for each interval for all subjects
createPredictData <- function(dataX, delta.upper){
    # Need to create a prediction set with the relevant data for this window
    IDInternal <- 1:nrow(dataX)
    n <- nrow(dataX)
    dataX$IDInternal <- IDInternal
    n.delta <- length(delta.upper)
    interval <- rep(1:n.delta, times=nrow(dataX))
    long.data <- dataX[rep(seq_len(nrow(dataX)), each=n.delta), ]
    long.data <- cbind(long.data, interval)
    # Sort the data by ID and interval
    long.data <- long.data[order(long.data$IDInternal, long.data$interval),]
    long.data$interval <- factor(long.data$interval)
    return(long.data)
}


# Function for computing predictions from a discrete time model --------
predBinSurv <- function(testDat, delta.upper, model, times) {
    #Format test data set for prediction into the same format as we used to fit the model
    testDat.pt <- createPredictData(dataX = testDat, delta.upper=delta.upper)

    #Get prediction for each of the models for conditional probabilities
    if(model$method %in% c("svmRadial", "cforest", "nnet")) {
        testDat.pt$pred <- predict(model, newdata=testDat.pt, type="prob")[,"yes"]
    } else {
        testDat.pt$pred <- predict(model, newdata=testDat.pt, type="prob")[,"1"]
    }

    #loop over individuals to calculate survival probabilities
    testSurvPreds = matrix(NA, nrow(testDat), length(times))
    rownames(testSurvPreds) <- 1:nrow(testDat)
    for(i in 1:length(times)) {
        delta.pred <- length(delta.upper[delta.upper<=times[i]]) #include all intervals where t_end <= t
        if(delta.pred==0) {
            testSurvPreds[,i] <- 0 #prob of experiencing the event at time 0 is 0
        } else {
            preds <- as.numeric(by(testDat.pt, testDat.pt[,"IDInternal"], function(x) {
                survP <- cumprod((1-x[["pred"]]))
                return(survP[delta.pred])
            }))
            testSurvPreds[,i] <- 1-preds
        }
    }
    return(testSurvPreds)
}


# Function for prediction from a logistic model ---------------------------
glmBinModel <- function(trainDat, timeVar = timeVar, statusVar = statusVar, bins=5, w, cens) {

    #create person-time data set
    ptData <- genPTdat(dat=trainDat, timeVar=timeVar, statusVar=statusVar, bins=bins, w=w, cens=cens)

    # creates a dataset with standard names, should use a unique name to avoid creating non unique names
    df.full <- subset(ptData, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))

    fitControl <- caret::trainControl(method="none", savePredictions="all")

    #Fit binary classification models
    modelGlm <- caret::train(N.delta ~ ., data=df.full, method="glm", na.action=na.omit,
                             trControl=fitControl, preProcess = c('center', 'scale'))

    return(list(delta.upper = unique(ptData[["delta.upper"]]),
                model = modelGlm))
}


# Function for prediction from a GAM --------------------------------------
gamBinModel <- function(trainDat, timeVar = timeVar, statusVar = statusVar, w, bins=5, cens) {

    #create person-time data set
    ptData <- genPTdat(dat=trainDat, timeVar=timeVar, statusVar=statusVar, bins=bins, w=w, cens=cens)

    # creates a dataset with standard names, should use a unique name to avoid creating non unique names
    df.full <- subset(ptData, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))

    fitControl <- caret::trainControl(method="none", savePredictions="all")

    #Fit binary classification models
    modelGam <- caret::train(N.delta ~ .,
                             data=df.full,
                             method="gam",
                             na.action=na.omit,
                             trControl=fitControl, preProcess = c('center', 'scale'),
                             tuneGrid = data.frame(method = "GCV.Cp", select = FALSE))

    return(list(delta.upper = unique(ptData[["delta.upper"]]),
                model = modelGam))
}


# Function for prediction from a glmnet -----------------------------------
glmnetBinModel <- function(trainDat, timeVar = timeVar, statusVar = statusVar, w, bins=5, cens,
                           alpha, lambda) {

    #create person-time data set
    ptData <- genPTdat(dat=trainDat, timeVar=timeVar, statusVar=statusVar, bins=bins, w=w, cens=cens)

    # creates a dataset with standard names, should use a unique name to avoid creating non unique names
    df.full <- subset(ptData, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))

    fitControl <- caret::trainControl(method="cv", number = 3, savePredictions="all", classProbs = F)

    #Fit binary classification models
    tuneGrid <- data.frame(alpha = alpha,
                           lambda = lambda)

    modelGlmnet <- caret::train(N.delta ~ .,
                                data=df.full,
                                method="glmnet",
                                na.action=na.omit,
                                trControl=fitControl,
                                tuneGrid = tuneGrid, preProcess = c('center', 'scale')
                                # tuneGrid = expand.grid(alpha = 1,lambda = seq(0.001, 0.1, by = 0.001))
    )

    return(list(delta.upper = unique(ptData[["delta.upper"]]),
                model = modelGlmnet))
}


# Function for prediction from a GBM --------------------------------------
gbmBinModel <- function(trainDat, timeVar = timeVar, statusVar = statusVar, w, bins=5, cens,
                        n.trees,
                        interaction.depth,
                        shrinkage,
                        n.minobsinnode) {

    #create person-time data set
    ptData <- genPTdat(dat=trainDat, timeVar=timeVar, statusVar=statusVar, bins=bins, w=w, cens=cens)

    # creates a dataset with standard names, should use a unique name to avoid creating non unique names
    df.full <- subset(ptData, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))

    tuneGrid <- data.frame(n.trees = n.trees,
                           interaction.depth = interaction.depth,
                           shrinkage = shrinkage,
                           n.minobsinnode = n.minobsinnode)

    fitControl <- caret::trainControl(method="none", savePredictions="all")

    modelGbm <- caret::train(N.delta ~ ., data=df.full, method="gbm", na.action=na.omit,
                             trControl=fitControl, verbose=FALSE, tuneGrid = tuneGrid,
                             preProcess = c('center', 'scale'))

    return(list(delta.upper = unique(ptData[["delta.upper"]]),
                model = modelGbm))
}


# Function for prediction from a SVM --------------------------------------
svmBinModel <- function(trainDat, timeVar = timeVar, statusVar = statusVar, w, bins=5, cens,
                        C, sigma) {

    #create person-time data set
    ptData <- genPTdat(dat=trainDat, timeVar=timeVar, statusVar=statusVar, bins=bins, w=w, cens=cens)

    # creates a dataset with standard names, should use a unique name to avoid creating non unique names
    df.full <- subset(ptData, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))

    fitControl <- caret::trainControl(method="none", savePredictions="all", classProbs =  TRUE)
    df.full$N.delta <- factor(df.full$N.delta, labels=c("no", "yes"))

    tuneGrid <- data.frame(C = C, sigma=sigma)

    modelSvm <- caret::train(N.delta ~ ., data=df.full, method="svmRadial",
                             na.action=na.omit, trControl=fitControl, verbose=FALSE,
                             tuneGrid = tuneGrid,
                             preProcess = c("center","scale"))

    return(list(delta.upper = unique(ptData[["delta.upper"]]),
                model = modelSvm))
}


# Function for prediction from a conditional inference forest -------------
cforestBinModel <- function(trainDat, timeVar, statusVar, w, bins=5, cens,
                            mtry) {

    #create person-time data set
    ptData <- genPTdat(dat=trainDat, timeVar=timeVar, statusVar=statusVar, bins=bins, w=w, cens=cens)

    # creates a dataset with standard names, should use a unique name to avoid creating non unique names
    df.full <- subset(ptData, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))

    fitControl <- caret::trainControl(method="none", savePredictions="all", classProbs =  TRUE)
    df.full$N.delta <- factor(df.full$N.delta, labels=c("no", "yes"))

    tuneGrid <- data.frame(mtry = mtry)

    modelCforest <- caret::train(N.delta ~ ., data=df.full, method="cforest", tuneGrid=tuneGrid,
                                 na.action=na.omit, trControl=fitControl, preProcess = c('center', 'scale'))

    return(list(delta.upper = unique(ptData[["delta.upper"]]),
                model = modelCforest))
}


# Function for prediction from a Neural network ---------------------------
nnetBinModel <- function(trainDat, timeVar, statusVar, w, bins=5, cens,
                         size, decay) {

    #create person-time data set
    ptData <- genPTdat(dat=trainDat, timeVar=timeVar, statusVar=statusVar, bins=bins, w=w, cens=cens)

    # creates a dataset with standard names, should use a unique name to avoid creating non unique names
    df.full <- subset(ptData, select= -c(IDInternal, timeInternal, eventInternal, delta.lower, delta.upper))
    df.full$N.delta <- factor(df.full$N.delta, labels=c("no", "yes"))

    fitControl <- caret::trainControl(method="none", savePredictions="all", classProbs =  TRUE)

    tuneGrid <- data.frame(size = size, decay = decay)

    modelNnet <- caret::train(N.delta ~ ., data=df.full, method="nnet", tuneGrid=tuneGrid,
                              na.action=na.omit, trControl=fitControl, trace = FALSE, preProcess = c('center', 'scale'))

    return(list(delta.upper = unique(ptData[["delta.upper"]]),
                model = modelNnet))
}


# Function for prediction from a Cox model --------------------------------
predictCox <- function(model, dat, w){
    sf <- survfit(model, newdata = dat)
    pred_surv <- summary(sf, times=w)$surv
    return(1-pred_surv)
}


# Function for prediction from an RSF model -------------------------------
predictRSF <- function (object, newdata, times, ...) {
    N <- NROW(newdata)
    class(object) <- c("rfsrc", "grow")
    S <- randomForestSRC::predict.rfsrc(object, newdata = newdata)$survival
    if(N == 1) S <- matrix(S, nrow = 1)
    Time <- object$time.interest
    p <- cbind(1, S)[, 1 + prodlim::sindex(Time, times),drop = FALSE]
    if(NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop("Prediction failed") #prediction fails if predictions cannot be made for every subject
    return(c(1 - p))
}
