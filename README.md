# autoSurv: An automated Machine Learning (AutoML) tool for survival analysis
Krithika Suresh (**maintainer**, <krithika.suresh@cuanschutz.edu>), Cameron Severn 

[[Paper]] | [[Tutorial]]

## Installation
The `autoSurv` functions are available for download from Github. 
```
library(devtools)
install_github("ksuresh17/autoSurv")
library(autoSurv)
```

## Introduction 
`autoSurv` can be used to train and evaluate the performance of survival prediction models applied to time-to-event data. Currently, the function implements the continuous-time models Cox proportional hazards and random survival forests, and various classification models (logistic regression, elastic net, support vector machines, gradient boosting machines, neural network) in a discrete-time framework. 

The existing functions perform the following steps: 

* Model training: The training/test data split can be specified within the function. The training data set is used to build the models and to tune the hyperparameters. Additionally, the user can specify whether cross-validation in the training set should be used for hyperparameter tuning. 
* Hyperparameter tuning: Machine learning algorithms often utilize hyperparameters in the model training process that must be tuned to identify the optimal values that will optimize predictive performance. Additionally, for the discrete-time models we tune the number of intervals used to create the person-time data set. Hyperparameter tuning is performed using Bayesian optimization. In the current function implementation, this step cannot be skipped with the aim that users can take advantage of the automated nature of this function to perform tuning. 
* Model testing: The models built on the training data (using the tuned hyperparameters parameters) are applied to the test data and time-dependent performance metrics for right-censored data are reported. Specifically, we consider a measure of discrimination, area under the curve (AUC), and an overall performance measure, Brier score. 
* New predictions: The trained models can then be used to obtain prediction for a new individual or set of individuals. Predictive performance can be computed for predictions on an independent data set. 

## Functions
There are currently two main functions implemented: `autoSurv()`, which is used to train and test discrete-time prediction models, and `predict_autoSurv()`, which is used to obtain predicted survival probabilities for new individuals, and can be used to assess predictive performance in an independent test set. 

## `autoSurv()` example: PBC data set

### Load and standardize data
The `pbc.dat` data set from the R package `survival` contains the survival time (labeled `survTime`), an event indicator (labeled `statusComp`), and 17 predictors. We standardize the covariates prior to building the prediction models.

```r 
data(pbc, package="survival")
pbc$statusComp <- ifelse(pbc$status == 2, 1, 0)
pbc$survTime <- pbc$time/365.25
pbc.dat <- pbc[,c(5,7:22)]
pbc.dat <- pbc.dat[complete.cases(pbc.dat),]

#Standardize covariates
pbc.dat_covs <- subset(pbc.dat, select = -c(statusComp,survTime))
pbc.dat_covs_process <- caret::preProcess(pbc.dat_covs, method=c("center", "scale"))
pbc.dat_covs <- predict(pbc.dat_covs_process, pbc.dat_covs)
pbc.dat <- cbind(pbc.dat[,c("survTime","statusComp")], pbc.dat_covs)
```

### `autoSurv()`: Fit survival models to PBC data

We use an 80/20 training/test split (`testProp=0.2`) and perform 5-fold (`cv=TRUE`, `cvFold=5`) cross-validation in the training data set to tune hyperparameters. 

We specify the prediction time of interest as 3 years (`times=3`). We consider a maximum of 25 intervals (`bins.upper=25`) for creating the discrete-time data set. The optimal number of intervals will be tuned during cross-validation using Bayesian optimization. We specify that censored observations should be treated as being censored in the same interval in which they were last observed (`cens="same"`).

The specification of the seed (`seed=1101`) allows for these results to be reproduced. 

The following models are fit (`trainModels`): 

* Continuous-time: Cox proportional hazards (`cox`) and Random survival forest (`rsf`)
* Discrete-time models: logistic regression (`glm`), elastic net (`glmnet`), gradient boosting machines (`gbm`), support vector machines (`svm`), and neural networks (`nnet`)

Some of these models require hyperparameters that are also tuned using Bayesian optimization. We use the default values for parameters related to the specifications of the Bayesian optimization search scheme (`init=10`, `iter=20`). Tuning is performed by minimizing the Brier score of the maximum of the prediction times specified (`times`) or by minimizing the integrated Brier score when the prediction times is specified as a vector of multiple times. 

Note that the following code will take a substantial amount of computation time (18 min on our system) due to the cross-validation and optimization of the hyperparameters. Of the available discrete-time models, `cforest` takes the longest to run and so has been excluded from the current list of models to reduce computation time. The optimization cannot be skipped in the current specification, but computation time can be reduced by reducing the number of times the Bayesian Optimization is repeated (`iter`). The printout indicates the steps in the optimization, which can be silenced (`verbose.opt=FALSE`).

```r
pbc.results <- autoSurv(timeVar = "survTime", 
                        statusVar = "statusComp", 
                        data = pbc.dat, 
                        times =  c(3),
                        trainModels = c("cox","rsf","glm","glmnet","gbm","svm","nnet"), #cforest
                        bins.upper = 25,
                        cens = "half",
                        testProp = 0.2, 
                        seed = 1101,
                        cv = TRUE,
                        cvFold = 5,
                        verbose.opt = FALSE
)
```

### View tuned hyperparameters for fitted models

We can view the tuned hyperparameters from the Bayesian optimization for each of the models. We can see the optimal number of intervals selected by tuning for each of the discrete-time models. 

```r 
pbc.results$tuned_params
```

### Assess predictive performance

For the tuned models, we examine the predictive performance in the test data set using the following time-dependent metrics that account for right-censoring:

* AUC (`$auc`): Ranges from 0 to 1, where values closer to 1 indicate better discrimination, and 0.5 indicates discriminative ability similar to chance 
* Brier Score (`Brier` column in `$brier`): Similar to mean squared error, where values closer to 0 indicate better overall predictive performance. The "Null model" indicates the performance of a model with no covariates and can be used as a benchmark for assessing the performance of the other models. 
* R-squared (`R2` column in `$brier`): Scaled Brier score computed as (Model Brier score)/(Null model Brier score), where values closer to 1 indicate better overall predictive performance. 
* Integrated Brier score (`IBS` column in `$brier`): Integrates the Brier score over the prediction times (`times`) specified in `autoSurv`. Will be 0 if there is only one prediction time specified. 

A data set of all performance metrics can also be obtained (`$metrics`). 

```r
pbc.results$auc

pbc.results$brier

pbc.results$metrics
```

### `predict_autoSurv()`: Obtain predictions for an individual 

We can obtain the predicted probabilities at a set of prediction times for a particular individual using the `predict_autoSurv()` function, and extract the predicted survival probabilities for any of the models that were fit. 

```r
times_pred <- seq(0,3,by=0.01)
pbc.testResults.id <- predict_autoSurv(pbc.results, 
                                    newdata=pbc.dat[4,], 
                                    times=times_pred, 
                                    timeVar="survTime", 
                                    statusVar="statusComp")
```





