#' Predict cell types by random forest models or logistic regression
#' 
#' @param train_data a cell x gene dataframe used for training, must have the same genes as test_data.
#' @param label.name the colname with labels in train_data
#' @param test_data a new cell x gene dataframe that user wants to predict cell types.
#' @param method Methods of prediction. Can be clm, rf_ranger, rf_rpart. Default is clm.
#' @return A object result from stat::predict
#' @export
predict_cell <- function(train_data, test_data, label.name, method = 'clm',
                         return.model = F) {
  stage <- train_data[, label.name]
  # scale data by features
  train_data <- train_data[,-which(colnames(train_data) %in% label.name)]
  if(!all.equal(colnames(test_data), colnames(train_data))) {
    print(dim(train_data))
    print(dim(test_data))
    stop('The genes from training data and test data are not equal.')
  }
  train_scaled <- as.data.frame(scale(train_data))
  train_scaled$Stage <- stage
  
  test_data <- as.data.frame(scale(test_data))
  
  if (method == 'clm') {
    model <- ordinal::clm(
      Stage ~ ., 
      data = train_scaled,
      link = "logit"  # 等同于logistic
    )
    summary(model)
    if (return.model) return(model)
    predictlist <- predict(clm_model, newdata = test_data, type = "class")
  
  }else if (method == 'rf_ranger') {
    model <- ranger::ranger(
      Stage ~ .,
      data = train_scaled,
      num.trees = 1000,
      mtry = floor(sqrt(ncol(train_scaled) - 1)),  # 特征抽样数
      splitrule = "gini",                          # 支持有序分裂
      importance = "impurity"
    )
    summary(model)
    if (return.model) return(model)
    predictlist <- predict(model, data = test_data, type = "response")
    
  }else if (method == 'rf_rpart') {
    model <- rpart(Stage ~ ., data = train_scaled, method = "class")
    summary(model)
    if (return.model) return(model)
    predictlist <- predict(model, newdata = test_data, type = "class")
  }else{stop('Argument method = *** is only provided: clm, rf_ranger or rf_rpart.')}
  return(predictlist)
}
