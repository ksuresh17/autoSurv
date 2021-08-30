#' autoSurv
#'
#' @param timeVar string corresponding to the variable name of the time-to-event
#' @param statusVar string corresponding to the variable name of the status indicator
#' @param idVar string corresponding to the variable name of the unique identifiers of each observation
#' @param data data frame containing covariates for prediction
#' @param times vector of time points at which probably of event is predicted
#' @param testProp proportion of observations to be randomly selected for the test set
#' @param seed integer random seed
#' @param init number of randomly chosen points to sample the target function before Bayesian Optimization fitting the Gaussian Process
#' @param iter total number of times the Bayesian Optimization is to be repeated
#'
#' @return
#' @export
#'
#' @examples
autoSurv <- function(timeVar,
                     statusVar,
                     idVar,
                     data,
                     times,
                     testProp = 0.2,
                     seed = 1,
                     init = 10,
                     iter = 20) {

    # split into train / test
    set.seed(seed)
    testIndexes <- sample(x = 1:nrow(data),size = nrow(data)*testProp)
    testDat <- data[testIndexes, ]
    trainDat <- data[-testIndexes, ]

    # create survival formula
    .form <- as.formula(paste("Surv(", timeVar,",", statusVar, ") ~ ."))

    # cox model
    coxphModel <- survival::coxph(.form,
                                  data=trainDat)

    preds <- list()
    for (w in times){
        pred_name <- paste("cox",w, sep = "_")
        preds[pred_name] <- list(cox_prob(coxphModel, dat=testDat, t0=0, w=w))
    }



    # Random Survival Forests
    rf_cv <- function(ntree,nodesize){
        rfsrcModel <- randomForestSRC::rfsrc(.form,
                                             ntree = ntree,
                                             nodesize = nodesize,
                                             data=trainDat)

        briers <- NULL
        for (w in times){
            pred.RSF <- predictRSF(rfsrcModel, newdata=testDat, times=w)

            temp_Score <- riskRegression::Score(list( "RSF" = pred.RSF),
                                                formula=as.formula(paste("Surv(", timeVar,",", statusVar, ") ~ 1")),
                                                data=testDat,
                                                times=w)
            briers <- c(briers,
                        temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "RSF"])
        }

        #neg_sum_briers <- -sum(briers)
        neg_int_briers <- -pracma::trapz(x = times, y = briers)

        list(Score = neg_int_briers,
             Pred = 0)
    }
    print("Optimizing RSF")

    OPT_Res <- rBayesianOptimization::BayesianOptimization(rf_cv,
                                                           bounds = list(ntree = c(300L, 1000L),
                                                                         nodesize = c(5L, 20L)),
                                                           init_grid_dt = NULL,
                                                           init_points = init,
                                                           n_iter = iter,
                                                           acq = "ucb",
                                                           kappa = 2.576,
                                                           eps = 0.0,
                                                           verbose = TRUE)

    rfsrcModel <- randomForestSRC::rfsrc(.form,
                                         ntree = OPT_Res$Best_Par[["ntree"]],
                                         nodesize = OPT_Res$Best_Par[["nodesize"]],
                                         data=trainDat)

    for (w in times){
        rf_name <- paste("rf",w, sep = "_")
        preds[rf_name] <- list(predictRSF(rfsrcModel, newdata=testDat, times=w))

    }

    # GLM discrete time model
    glm_cv <- function(intervals){

        briers <- NULL
        for (w in times){
            pred.GLM <- BinModel(testDat = testDat, trainDat = trainDat,
                                 timeVar = "survTime", eventVar="statusComp", excludeVars = c("id","status","fold"),
                                 w=w, numInt=intervals, methods.lib="glm")


            temp_Score <- riskRegression::Score(list( "GLM" = pred.GLM),
                                                formula=as.formula(paste("Surv(", timeVar,",", statusVar, ") ~ 1")),
                                                data=testDat,
                                                times=w)
            briers <- c(briers,
                        temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "GLM"])
        }

        #neg_sum_briers <- -sum(briers)
        neg_int_briers <- -pracma::trapz(x = times, y = briers)

        list(Score = neg_int_briers,
             Pred = 0)
    }

    print("Optimizing DiscreteTime-GLM")

    OPT_Res_GLM <- rBayesianOptimization::BayesianOptimization(glm_cv,
                                                              bounds = list(intervals = c(2L, 10L)),
                                                              init_grid_dt = NULL,
                                                              init_points = init,
                                                              n_iter = iter,
                                                              acq = "ucb",
                                                              kappa = 2.576,
                                                              eps = 0.0,
                                                              verbose = TRUE)


    for (w in times){
        glm_name <- paste("glm",w, sep = "_")
        preds_GLM <- c(BinModel(testDat = testDat,
                                trainDat = trainDat,
                                timeVar = timeVar,
                                eventVar = statusVar,
                                excludeVars = c("id","status","fold"),
                                w=w,
                                numInt=OPT_Res_GLM$Best_Par[["intervals"]],
                                methods.lib="gbm"))
        preds[glm_name] <- list(preds_GLM)

    }

    # GBM discrete time model
    gbm_cv <- function(intervals){

        briers <- NULL
        for (w in times){
            pred.GBM <- BinModel(testDat = testDat, trainDat = trainDat,
                                 timeVar = "survTime", eventVar="statusComp", excludeVars = c("id","status","fold"),
                                 w=w, numInt=intervals, methods.lib="gbm")


            temp_Score <- riskRegression::Score(list( "GBM" = pred.GBM),
                                                formula=as.formula(paste("Surv(", timeVar,",", statusVar, ") ~ 1")),
                                                data=testDat,
                                                times=w)
            briers <- c(briers,
                        temp_Score$Brier$score$Brier[temp_Score$Brier$score$model == "GBM"])
        }

        #neg_sum_briers <- -sum(briers)
        neg_int_briers <- -pracma::trapz(x = times, y = briers)

        list(Score = neg_int_briers,
             Pred = 0)
    }

    print("Optimizing DiscreteTime-GBM")

    OPT_Res_GBM <- rBayesianOptimization::BayesianOptimization(gbm_cv,
                                                               bounds = list(intervals = c(2L, 10L)),
                                                               init_grid_dt = NULL,
                                                               init_points = init,
                                                               n_iter = iter,
                                                               acq = "ucb",
                                                               kappa = 2.576,
                                                               eps = 0.0,
                                                               verbose = TRUE)


    for (w in times){
        gbm_name <- paste("gbm",w, sep = "_")
        preds_GBM <- c(BinModel(testDat = testDat,
                                trainDat = trainDat,
                                timeVar = timeVar,
                                eventVar = statusVar,
                                excludeVars = c("id","status","fold"),
                                w=w,
                                numInt=OPT_Res_GBM$Best_Par[["intervals"]],
                                methods.lib="gbm"))
        preds[gbm_name] <- list(preds_GBM)

    }



    temp_Score <- riskRegression::Score(preds,
                                        formula=as.formula(paste("Surv(", timeVar,",", statusVar, ") ~ 1")),
                                        data=testDat,
                                        times=times)

    scoreAUC <- temp_Score$AUC$score
    scoreAUC <- scoreAUC[stringr::str_extract(scoreAUC$model,"[0-9]") == scoreAUC$times,]
    scoreBrier <- temp_Score$Brier$score
    scoreBrier <- scoreBrier[stringr::str_extract(scoreBrier$model,"[0-9]") == scoreBrier$times,]
    scoreBrier <- scoreBrier[!is.na(scoreBrier$model),]

    list("auc" = scoreAUC, "brier" = scoreBrier, "pred_probabilities" = preds,"models" = list("bestRF" = rfsrcModel, "bestCox" = coxphModel))

}

createDiscreteIntervals2 <- function(time, event, numInt) {
    n.delta <- min(numInt, length(unique(time[event==1])))
    probs.delta <- seq(from=0, to=1, length.out=(n.delta+1))
    delta.upper <- quantile(time[event==1], probs=probs.delta, names=FALSE)
    return(delta.upper[-1])
}

## This function was adapted from code from the Github account of Eric Polley
## https://github.com/ecpolley/SuperLearner_Old/blob/master/R/createDiscrete.R
## Also see Chapter 16 "Super Learning for Right-Censored Data" in the book
## "Targeted Learning: Causal Inference for Observational and Experimental Data"
## by van der Laan and Rose, Springer, 2011.

# Function that returns binary process data set ---------------------------
createDiscreteDat <- function(time, event, dataX, delta.upper) {
    n.time <- length(time)
    #Assumes t0=0
    delta.lower <- c(0, delta.upper[-length(delta.upper)])
    n.delta <- length(delta.upper)
    ID <- 1:nrow(dataX)
    dat_i <- cbind(ID, time, event, dataX)
    interval <- rep(1:length(delta.upper), times=nrow(dataX))

    #Create a long data set that repeats each persons data n.delta (number of intervals) times
    long.dat <- dat_i[rep(seq_len(nrow(dat_i)), n.delta), ]

    N.delta <- rep(NA, nrow(long.dat))
    long.dat <- cbind(long.dat, delta.lower, delta.upper, N.delta, interval)

    # Include censored people in the interval in which they were censored
    long.dat$N.delta <- ifelse(long.dat$time > long.dat$delta.upper, 0,
                               ifelse(long.dat$event==1, ifelse(long.dat$time <= long.dat$delta.lower, NA, 1),
                                      ifelse(long.dat$time>long.dat$delta.lower, 0, NA)))

    m <- delta.upper[n.delta]
    long.dat$N.delta <- ifelse(long.dat$time == m & long.dat$delta.upper == m,
                               ifelse(is.na(long.dat$N.delta), 0, long.dat$N.delta), long.dat$N.delta)

    # Drops the intervals to at which the individual no longer contributes
    long.dat <- long.dat[!is.na(long.dat$N.delta), ]
    # Sort the data
    long.dat <- long.dat[order(long.dat$ID, long.dat$delta.lower),]
    # Set interval as a factor
    long.dat$interval <- as.factor(long.dat$interval)
    # Set outcome as a factor
    long.dat$N.delta <- as.factor(long.dat$N.delta)
    return(long.dat)
}

# Function to create person-time data set ---------------------------------
# dat: (data.frame) data set to convert (with columns named "survTime", "event")
# numInt: (integer) number of discrete intervals (hyperparameter that can be tuned)
# w: (integer) prediction interval
genPTdat <- function(dat, timeVar = "survTime", eventVar = "event", numInt = 5, w) {
    time.ind <- which(eval(timeVar) == colnames(dat))
    event.ind <- which(eval(eventVar) == colnames(dat))

    dat_i <- dat
    dat_i[dat_i$survTime > w, event.ind] <- 0 #ignore events after the prediction window
    dat_i[dat_i$survTime > w, time.ind] <- w #administratively censor at w

    # Identify the times t1, t2,... tJ (assumes t0=0)
    delta.upper <- createDiscreteIntervals2(time = dat_i[,time.ind],
                                            event = dat_i[event.ind],
                                            numInt = numInt)

    # Creates a new data set where each person has a row corresponding to the discrete time intervals
    # 0: in intervals where alive but does not have event
    # 1: in intervals where experience the event
    # delta.upper and delta.lower columns: are the same for everyone and are the discretized time intervals based on quantiles
    dat_i.X <- createDiscreteDat(time = dat_i[,time.ind],
                                 event = dat_i[,event.ind],
                                 dataX = dat_i[,-c(time.ind,event.ind)],
                                 delta.upper = delta.upper)

    return(dat_i.X)
}

# Function that formats data we will make predictions on ------------------
#Similar to above function but does not create outcome variable and creates n.delta rows for all subjects
createPredictData2 <- function(dataX, time){
    # Need to create a prediction set with the relevant data for this window
    n <- nrow(dataX)
    n.delta <- length(time)
    interval <- rep(1:n.delta, times=nrow(dataX))

    long.data <- dataX[rep(seq_len(nrow(dataX)), each=n.delta), ]
    long.data <- cbind(long.data, interval)
    long.data$interval <- factor(long.data$interval)
    long.data <- as.data.frame(long.data)
    return(long.data)
}


BinModel <- function(testDat, trainDat, timeVar = "survTime", eventVar = "statusComp", excludeVars = NULL, w, numInt=5, methods.lib=c("gbm","glm")) {

    if(!is.null(excludeVars)) {
        testDat2 <- testDat[,-which(colnames(testDat) %in% eval(c("id","status","fold")))]
        trainDat2 <- trainDat[,-which(colnames(trainDat) %in% eval(c("id","status","fold")))]
    }

    #create person-time data set
    ptData <- genPTdat(dat=trainDat2, timeVar=timeVar, eventVar=eventVar, numInt=numInt, w=w)

    df.X <- subset(ptData, select= -c(ID, time, event, delta.lower, delta.upper, N.delta))
    df.full <- subset(ptData, select= -c(ID, time, event, delta.lower, delta.upper))

    fitControl <- caret::trainControl(method="none", savePredictions="all")

    #Fit binary classification models
    if("gbm"%in%methods.lib){
        modelGbm <- caret::train(N.delta ~ ., data=df.full, method="gbm", na.action=na.omit, trControl=fitControl, verbose=FALSE)
    }
    if ("glm"%in%methods.lib) {
        modelGlm <- caret::train(N.delta ~ ., data=df.full, method="glm", na.action=na.omit, trControl=fitControl)
    }

    ids.test <- 1:nrow(testDat2)
    delta.bound <- sort(unique(ptData["delta.lower"][ptData["delta.lower"]!=0] ))
    nhor <- length(delta.bound)

    ## Must put the data for prediction into the same format as we used to fit the model
    testDat.pt <- createPredictData2(dataX = testDat2, time=delta.bound)

    #Get prediction for each of the models for conditional probabilities
    testDat.pred <- testDat.pt[1]
    if("gbm"%in%methods.lib){
        testDat.pred$gbm <- predict(modelGbm, newdata=testDat.pt, type="prob")[,"1"]
    }
    if("glm"%in%methods.lib){
        testDat.pred$glm <- predict(modelGlm, newdata=testDat.pt, type="prob")[,"1"]
    }

    #loop over individuals to calculate survival probabilities
    testSurvPreds = matrix(NA, length(ids.test), length(methods.lib))
    rownames(testSurvPreds) <- ids.test
    #better way to do this by id instead of relying on data set being correctly formulated?
    for(j in 1:length(ids.test)){
        beg <- ((j-1)*nhor) +1
        en <- j*nhor
        for(k in 1:(ncol(testDat.pred)-1)) {
            survP <- cumprod((1-testDat.pred[beg:en, k+1]))
            testSurvPreds[ids.test[j],k] <- survP[nhor]
        }
    }
    return(1-testSurvPreds)
}

# Function for prediction from a Cox model --------------------------------
cox_prob <- function(model, dat, t0, w){
    ret <- NULL
    for(i in 1:nrow(dat)) {
        sf <- survfit(model, newdata = dat[i,]) #get the predicted cummulative hazard for each individual based on fit model
        sf2 <- data.frame(time = sf$time, surv = sf$surv, Haz = -log(sf$surv))
        tmp <- dynpred::evalstep(sf2$time,sf2$Haz,c(t0,t0+w),subst=0) #landmark time and event horizon defined here to evaluate CHF
        Fw <- exp(-(tmp[2]-tmp[1])) #this is calculating probability of event?
        ret <- c(ret, Fw)
    }
    return(1-ret)
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
