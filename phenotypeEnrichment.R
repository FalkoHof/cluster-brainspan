
calculateEnrichment <- function(kMeansClusters, tfData, idPhenotypes, filename='outputData/enchrichment.csv'){
  #create a data frame where the enrichment results will be stored in 
  resultsData <- data.frame(stringsAsFactors=F)

  clusterData <- data.frame(ensembl_gene_id=names(kMeansClusters$cluster),cluster_id=kMeansClusters$cluster, stringsAsFactors=F)
  
  #filter the frame for transciption factors
  clusterTFsData <- clusterData[c(tfData$Ensembl.ID),]
  clusterTFsData <- clusterTFsData[complete.cases(clusterTFsData),]

  #save the data
  write.csv(clusterData,  file='outputData/clusters.csv') 
  write.csv(clusterTFsData,  file='outputData/tfClusters.csv') 
  
  #create a frame containg the clusters with >=1 TF and merge it with the phenotype data
  clustersOfInterest <- merge(clusterData, phenotypeData, by.x="ensembl_gene_id", by.y="ensembl_gene_id")
  clustersOfInterest <- unique(clustersOfInterest)
  #get the cluster id of the clusters containing >=1 TF
  indexesClustersOfInterest <-unique(clusterTFsData$cluster_id)
  
  print(paste('Clusters of Interest:', indexesClustersOfInterest))
  
  #loop over the clusters with >=1 TF
  for (x in indexesClustersOfInterest){
    #subset the clusters of interest into  
    oneCluster<- clustersOfInterest[clustersOfInterest$cluster_id == x,]
    #get phenotypes for the TFs
    phenotypesTFs <-oneCluster[oneCluster$is_TF,c("sysid_additional_class_type","ensembl_gene_id","cluster_id")]
    phenotypesHypergeometric <-unique(phenotypesTFs$sysid_additional_class_type)

    #reduce data frame to gene name and phenotype --> kick out potential multiple ids mapping to the same gene 
    oneCluster$ensembl_gene_id <- NULL
    oneCluster$cluster_id <- NULL
    oneCluster <- unique(oneCluster)
        
    #calculate the phenotype enrichment via a hypergeometric distribution
    #q, number of white balls drawn without replacement from an urn which contains both black and white balls.
    #m, the number of white balls in the urn.
    #n, the number of black balls in the urn.
    #k, the number of balls drawn from the urn.
  
    for (p in phenotypesHypergeometric){
      # get number genes within one cluster that have phenotype p
      q <- dim(oneCluster[oneCluster$sysid_additional_class_type == p,])[1]
      # get total number of genes that have phenotype p
      m <- dim(idPhenotypes[idPhenotypes$sysid_additional_class_type == p,])[1]
      # get tota number of genes that dont have phenotype p
      n <- dim(idPhenotypes[idPhenotypes$sysid_additional_class_type!=p,])[1]
      # get the number of phenotypes known for the cluster mebers
      k <- dim(oneCluster)[1]
    
      pValue <- phyper(q,m,n,k, lower.tail=FALSE)
      
      if(q>k)
        warning(paste('q>k, invalid parameters', 'q:', q, 'k:', k, 'phenotype:', p))
      if(q>m)
        warning(paste('q>m, invalid parameters', 'q:', q, 'k:', k, 'phenotype:', p))  
  
      resultsData[x,p] <- pValue
      #resultsData$pValues[x,p] <- pValue
    }
  }
  
  #bonfPValue <- 0.01/sum(!is.na(resultsData$pValues))
  #resultsData$pValues_bonferroni <- resultsData$pValues[resultsData$pValues <=bonfPValue]
  
  write.csv(resultsData,  file=filename)
  return(resultsData)
}