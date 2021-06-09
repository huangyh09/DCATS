
#' Calculate stochastic transition matrix between clusters from a data frame 
#' including information about clustering
#'
#' The transition probability from cluster i to j is calculated based on the 
#' information used to cluster cells. It is estimated by the misclassification 
#' rate from cluster i to j comparing the original labels with the labels predicted 
#' by support vector machine with 5-fold cross validation.
#'
#' @param dataframe a data frame contains the information used for clustering 
#' and the original label of each cell. The original lable should have the 
#' column name `clusterRes`.
#' 
#' @export
#'
#' @note `tidyverse` and `tidymodels` need to be load to allow this function to work

svm_simMat = function(dataframe){
    cv <- vfold_cv(dataframe, v = 5)
    recipe <- recipe(clusterRes ~ ., data = dataframe)
      model <-
        svm_rbf() %>% 
        set_mode("classification") %>%
        set_engine("kernlab")
    workflow <- workflow() %>%
      add_recipe(recipe) %>%
      add_model(model)
    predDF = data.frame()
    for (i in 1:5) {
      onefold_split = cv[[1]][[i]]
      fit <- workflow %>%
        last_fit(onefold_split)
      pred = fit %>% 
        collect_predictions() %>% 
        mutate(pred = .pred_class) %>% 
        select(pred, clusterRes)
      predDF = rbind(predDF, pred)
    }
    conf.mat <- table(predDF$clusterRes, predDF$pred)
    simil_mat <- t(t(conf.mat)/apply(conf.mat,2,sum))
  return(simil_mat)
}
