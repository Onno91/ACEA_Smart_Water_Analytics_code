model_handmade_folds <- function(data, horizon, dataset, lag, features, basefeatures){
  
  #basefeatures <- 'Depth'
  
  # Make lags:
  features <- grep(features,names(data),value = T)
  basefeatures <- grep(basefeatures,names(data),value = T)
  
  for(i in 1:length(features)){
    for(j in 1:lag){
      data$temp <- Lag(data[,features[i],+j])
      names(data)[which(names(data)=='temp')] <- paste(features[i],j, sep = '_')
    }
  }
  
  
  data <- data[,which(colMeans(!is.na(data))>.2)]
  
  # Inlude seasonality:
  data$year <- as.numeric(substr(data$Date,7,10))
  data$year <- data$year - min(data$year) + 1
  data$month <- as.numeric(substr(data$Date,4,5))
  data$quarter <- ifelse(data$month <= 3,1,
                         ifelse(data$month >=4 & data$month <= 6,2,
                                ifelse(data$month >=7 & data$month <= 9,3,
                                       ifelse(data$month >9,4,NA))))
  
  data_long <- tidyr::pivot_longer(data, basefeatures,names_to = 'location', values_to = 'depth_to_groundwater')
  
  data_long <- data_long[complete.cases(data_long),]
  data_long <- data_long[which(data_long$depth_to_groundwater != 0),]
  
  #data_model <- data_long[,-grep('location|Date|name',names(data_long))]
  
  temp <- data_long[,which(!names(data_long)%in%c('depth_to_groundwater','Date','location'))]
  nzv <- nearZeroVar(temp)                                                       # excluding variables with very low frequencies
  if(length(nzv)>0){temp <- temp[, -nzv]}
  i <- findCorrelation(cor(temp))                                                # excluding variables that are highly correlated with others
  if(length(i) > 0) temp <- temp[, -i]
  i <- findLinearCombos(temp)                                                    # excluding variables that are a linear combination of others
  if(!is.null(i$remove)) temp <- temp[, -i$remove]
  data_model <- data_long[,c('depth_to_groundwater','Date','location', names(temp))]
  
  data_model$Date <- as.Date(as.character(data_model$Date), format = '%d/%m/%Y')
  
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
  
  for(i in 1:length(locations)){
    for(j in 1:length(output_final)){
      if(length(locations)==1){
        
        output_temp <- index_hand_design(data_model,output_final[[j]], locations[i], horizon, location_one = 'yes') 
      } else {
        output_temp <- index_hand_design(data_model,output_final[[j]], locations[i], horizon)
      }
      if(j == 1){
        final_inner <- output_temp
      } else {
        final_inner <- c(final_inner, output_temp)
      }
    }
    if(i == 1){
      final <- final_inner
    } else {
      final <- c(final, final_inner)
    }
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
