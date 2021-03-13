LLTO_Model <- function(data, name_of_dataset, lag, time_out_horizon, features, predictors){
  ## Function purpose: Use Leave Location and Time Out Moddeling (LLTO) on the 
  ## DataSets of the Kaggle ACEA Smart Water Analytics Challenge. Datasets are 
  ## time series with multiple predictors and features. Lags are used as features.
  ## A horizon gets specified to model the Time Out Cross validation.
  
  ## Function Output: Train Object
  
  ## Function input:
  # data                Dataset from Kaggle Acea Smart Water Analytics challenge, 
  #                     containing Date, numeric predictors & feature set
  # name_of_dataset     Name of input dataset, used to save trained model
  # lag                 Amount of lags used for features (integer input)
  # time_out_horizon    Horizon of folds used in Time Out Cross validation
  # features            Names of features used in dataset (numeric variables)
  # predictors          Names of predictors used in dataset (numeric variables)
  
  ## Variable formatting:
  # Date variable:
  data_model$Date <- as.Date(as.character(data_model$Date), format = '%d/%m/%Y')
  
  ## Feature addition:
  # Creating lags for i features in j lags
  for(i in 1:length(features)){
    for(j in 1:lag){
      data$temp <- quantmod::Lag(data[,features[i],+j])
      names(data)[which(names(data)=='temp')] <- paste(features[i],j, sep = '_')
    }
  }
  # Inlude seasonality:
  data$year    <- as.numeric(substr(data$Date,7,10))
  data$year    <- data$year - min(data$year) + 1
  data$month   <- as.numeric(substr(data$Date,4,5))
  data$quarter <- ifelse(data$month <= 3,1,
                         ifelse(data$month >=4 & data$month <= 6,2,
                                ifelse(data$month >=7 & data$month <= 9,3,
                                       ifelse(data$month >9,4,NA))))  
  # remove uncommon newly created features:
  data <- data[,which(colMeans(!is.na(data))>.2)]
  
  ## Data shaping for modelling
  # Long format by predictors for LLTO modelling:
  data_long <- tidyr::pivot_longer(data, predictors,names_to = 'location', values_to = 'depth_to_groundwater')
  
  ## Remove missing variables:
  # Complete cases
  data_long <- data_long[complete.cases(data_long),]
  
  ## Statistical preprocessing:
  # Remove outliers
  data_long <- data_long[which(data_long$depth_to_groundwater != 0),]
  # Preprocess features for modelling
  temp <- data_long[,which(!names(data_long)%in%c('depth_to_groundwater','Date','location'))]
  # excluding variables with very low frequencies
  nzv <- nearZeroVar(temp)                                                       
  if(length(nzv)>0){temp <- temp[, -nzv]}
  # excluding variables that are highly correlated with others
  i <- findCorrelation(cor(temp))                                               
  if(length(i) > 0) temp <- temp[, -i]
  # excluding variables that are a linear combination of others
  i <- findLinearCombos(temp)                                                   
  if(!is.null(i$remove)) temp <- temp[, -i$remove]
  # Selection of final feature set
  data_model <- data_long[,c('depth_to_groundwater','Date','location', names(temp))]
  
  ## 
  
  # Handmade indexes:
  index_hand_design <- function(data,period, location, horizon, location_one = NULL){
    horizon2 <- max(period)-horizon
    if(!is.null(location_one)){
      indexin <- which(data$Date >= min(period) & data$Date <= horizon2)
      indexout <- which(data$Date > horizon2 & data$Date <= max(period))
      
    } else {
      indexin <- which(data$Date >= min(period) & data$Date <= horizon2 & data$location != location)
      indexout <- which(data$Date > horizon2 & data$Date <= max(period) & data$location == location)
    }
    output <-c(list(indexin),list(indexout))
    output
  }
  
  periods <- round(length(seq.Date(from = min(data_model$Date),to = max(data_model$Date), by = 'day'))/horizon,0)
  dates   <- seq.Date(from = min(data_model$Date),to = max(data_model$Date), by = 'day')
  indices <- 1:periods*horizon
  
  periods_final <- dates[indices]
  periods_final <- periods_final[!is.na(periods_final)]
  
  stopifnot(length(periods_final)>=4)
  
  for(i in 3:length(periods_final)){
    output <- list(c(periods_final[i-2], periods_final[i]))
    if(i <= 3){
      output_final <- output
    } else {
      output_final <- c(output_final, output)
    }
  }
  
  locations <- unique(data_model$location)
  
  if(length(locations)==1){
    final_inner <- lapply(1:length(output_final), function(x){index_hand_design(data_model,
                                                                                output_final[[x]],
                                                                                locations, 
                                                                                horizon,
                                                                                location_one = 'yes')})
  }
  if(length(locations)>1){
    # Create dataframe with possible combinations:
    sequences <- data.frame(location = rep(locations,length(output_final)), 
                                           lengths = rep(1:length(output_final), length(locations)))
    final <- mapply(index_hand_design, data = data_model, sequences[,2], sequences[,1], horizon = horizon)
  }
    
  
  
  index_final <- list(index = final[seq(1, length.out = length(locations)*length(output_final), by = 2)], 
                      indexOut = final[seq(2, length.out =length(locations)*length(output_final), by = 2)])
  
  
  fitcontrol <- trainControl(verboseIter = T,
                             index = index_final$index,
                             indexOut = index_final$indexOut)
  
  gbmGrid <-  expand.grid(interaction.depth = c(1,2,4), 
                          n.trees = 1:4000, 
                          shrinkage = c(0.01), 
                          n.minobsinnode = c(2,5))
  
  if(length(locations)>1){
    for(i in 1:length(locations)){
      data_model$temp <- ifelse(data_model$location == locations[i],1,0)
      names(data_model)[which(names(data_model) == 'temp')] <- paste(locations[i],'ind', sep = '_')
    }
  }
  
  err <- try(load(paste(maindir, modeldir, paste('data_model = ',dataset,'horizon =',horizon,'.RData', sep = ''), sep = '/')))
  if(err != 'train'){
    
    train <- train(depth_to_groundwater ~ . , method = 'gbm', trControl = fitcontrol, tuneGrid = gbmGrid, 
                   data = data_model[,-grep('Date|location',names(data_model))])
    save(train, file = paste(maindir, modeldir, paste('data_model = ', dataset,'horizon =', horizon,'.RData', sep = ''), sep = '/'))
  }
  return(train)
}
