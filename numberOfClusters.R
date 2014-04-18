###########################################
# determine number of clusters - gap method
###########################################

gapMethod <- function(maxNClusters, data, iterMax = 30, nStart = 100, b = 100){
  require(cluster)
  print('Starting Gap method...')
  return(clusGap(data, kmeans, maxNClusters, iter.max = iterMax, nstart = nStart, B = b, verbose = interactive()))
}





#############################################################
# determine number of clusters - within groups sum of squares
#############################################################

wssMethod <- function(maxNClusters, data, filenamePlot ='outputData/K-meansSSE.png', filenameRaw='outputData/K-meansSSE.png', writeRaw = FALSE, numberOfCPUCores=3){
  require(doMC)
  registerDoMC(numberOfCPUCores)
  
  print('WSS method... - Starting')
  
  #caclulate and plot wss as second controll for the number of clusters
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  wssPara <- foreach(i = 2:maxNClusters) %dopar%{sum(kmeans(data,centers=i, iter.max=30, nstart=100)$withinss)}
  wss <- c(wss,wssPara)
  
  #plot wss to file 
  png(file=filenamePlot)
  plot(1:length(wss), wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  dev.off()
  
  #write raw data to file
  if(writeRaw)
    lapply(wss, write, file=filenameRaw, append=TRUE)
  
  print('WSS method... - Done')
  return(wss)
}


