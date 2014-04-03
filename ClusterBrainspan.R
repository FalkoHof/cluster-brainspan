library('biomaRt')
library(mclust)
library(doMC)
library(cluster)

#register 3 cores for parallelization
registerDoMC(3)


########################
# load and process data:
########################

#load data
print('loading required data...')
data <- read.csv('/inputData/brainspan/expression_matrix.csv', header=F,stringsAsFactors=F)[,(-1)]
column_metadata<- read.csv('/inputData/brainspan/columns_metadata.csv')
rows_metadata <- read.csv('/inputData/brainspan/rows_metadata.csv')
tfData <- read.csv('/inputData/TF_List_Nature.txt', header=T, stringsAsFactors=F, sep='\t')
phenotypeData <-read.csv('/inputData/databaseAnnettePhenotypes.csv',header=T, stringsAsFactors=F, fill=T, strip.white=TRUE)
print('loading required data... - Done')

#add lables
print('labeling brainspan data...')
rownames(data) <- rows_metadata$ensembl_gene_id
colnames(data) <- paste(column_metadata$structure_acronym, column_metadata$age, column_metadata$donor_name, sep="-")
print('labeling brainspan data... - Done')

#kick old, partly missing ensembl ids and mitochondrial genes
print('formatting phenotype data...')
phenotypeData <- phenotypeData[phenotypeData$inheritance_pattern !='Mitochondrial',]
phenotypeData$ensembl_id <- NULL
phenotypeData$hgnc_id <- NULL
phenotypeData$gene_groups<- NULL
phenotypeData$inheritance_pattern <- NULL

#map entrez ids+phenotypes to ensembl ids
print('Mapping data via biomart...')
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
mapping <- getBM(attributes=c('ensembl_gene_id','entrezgene'), filters = "entrezgene", values = phenotypeData$entrez_id, mart = ensembl, uniqueRows=F)
print('Mapping data via biomart... - Done')

#merge the mapping and the phenotype data
phenotypeData <- merge(phenotypeData, mapping, by.x="entrez_id", by.y= "entrezgene")
write.csv(phenotypeData,  file='/outputData/phenotypeDataWEnsemble.csv') 

#select brainspan data of id genes and remove incomplete data
print('Subsetting brainspan data to ID genes and adding phenotype...')
subsetData <- data[c(mapping$ensembl_gene_id),]
subsetData <- subsetData[complete.cases(subsetData),]
write.csv(subsetData,  file='/outputData/brainspanIDGenes.csv') 

#kick out unneeded data from the phenotype data frame size can be used as indicator for the number of present phenotype
phenotypeData$omim_gene_id <- NULL
phenotypeData$omim_gene_id <- NULL
phenotypeData$hprd_id <- NULL
phenotypeData$inheritance_type<- NULL
phenotypeData$entrez_id <-NULL
phenotypeData$sysid_additional_class_type <-NULL

phenotypeData$is_TF <- phenotypeData$ensembl_gene_id %in% tfData$Ensembl.ID
phenotypeData <- phenotypeData[complete.cases(phenotypeData),]

totalPhenotypes <- unique(data.frame(official_gene_symbol = phenotypeData$official_gene_symbol, sysid_main_class_type= phenotypeData$sysid_main_class_type))

###########################################
# clustering -BIC of gaussian mixed models
###########################################
print('Start clustering...') 
#do the clustering 
cluster<-Mclust(subsetData)
nOfClusters <- dim(cluster$z)[2]
#write some output files
write.csv(cluster$z,  file='/outputData/mclustProbablities.csv') 
write.csv(cluster$classification,  file='/outputData/mclustClassification.csv') 
print('Start clustering... - Done') 


####################################################Start alternativ methods for clustering and determining k #########################################
####################################################Can be skipped#####################################################################################
###########################
# clustering - gap method
###########################
maxNClusters <- 30
myKmeans <- function(x, centers) {
  print(paste('Percent complete: ', round(((centers/maxNClusters)*100), digits=2), '%'))
  return(kmeans(x, centers, iter.max = 30, nstart = 10)) 
}
gapMethod <- clusGap(subsetData, myKmeans, 30, B = 100, verbose = interactive())


############################################
# clustering - within groups sum of squares
############################################

#caclulate and plot wss as second controll for the number of clusters
wss <- (nrow(subsetData)-1)*sum(apply(subsetData,2,var))
wssPara <- foreach(i = 2:maxNClusters) %dopar%{sum(kmeans(subsetData,centers=i, iter.max=30, nstart=10)$withinss)}
wss <- c(wss,wssPara)

png(file='/outputData/K-meansSSE.png')
plot(1:length(wss), wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

lapply(wss, write, file='/outputData/K-meansSSEValues.txt', append=TRUE) 

print('Dermining number of clusters per sample... - Done')



##############################################
# clustering - k means with calculated value
##############################################

#do the clustering 
print('Performing k-means clustering...')
kMeansClusters <- kmeans(subsetData ,nOfClusters, iter.max = 30, nstart = 1000)

print('Performing k-means clustering... Done')


####################################################End alternativ methods for clustering and determining k #########################################



#################################
# test for phenotype enrichment:
#################################

#get the clustering results
clusterData <- data.frame(ensembl_gene_id=names(cluster$classification),cluster_id=cluster$classification, stringsAsFactors=F)

#filter the frame for transciption factors
tfEnsembl <- tfData$Ensembl.ID
clusterTFsData <- clusterData[c(tfEnsembl),]
clusterTFsData <- clusterTFsData[complete.cases(clusterTFsData),]
#save the data
write.csv(clusterData,  file='/outputData/clusters.csv') 
write.csv(clusterTFsData,  file='/outputData/tfClusters.csv') 

#create a frame containg the clusters with >=1 TF and merge it with the phenotype data
clustersOfInterest <- merge(clusterData, phenotypeData, by.x="ensembl_gene_id", by.y= "ensembl_gene_id")
clustersOfInterest <- unique(clustersOfInterest)
#get the cluster id of the clusters containing >=1 TF
indexesClustersOfInterest <-unique(clusterTFsData$clusterID)

#create a data frame where the enrichment results will be stored in 
resultsData <- data.frame(row.names=clusterIndexes, stringsAsFactors=F)

#loop over the clusters with >=1 TF
for (x in clusterIndexes){
  #subset the clusters of interest into  
  oneCluster<- clustersOfInterest[clustersOfInterest$cluster_id == x,]
  #get phenotypes for the TFs
  phenotypesTFs <-oneCluster[oneCluster$is_TF,c("sysid_main_class_type","ensembl_gene_id","cluster_id")]
  phenotypesTFs <-unique(phenotypesTFs)
  phenotypesHypergeometric <-unique(phenotypesTFs$sysid_main_class_type)
  
  #calculate the phenotype enrichment via a hypergeometric distribution
  #q, number of white balls drawn without replacement from an urn which contains both black and white balls.
  #m, the number of white balls in the urn.
  #n, the number of black balls in the urn.
  #k, the number of balls drawn from the urn.
  
  for (p in phenotypesHypergeometric){
    # get number genes within one cluster that have phenotype p
    q <- dim(oneCluster[oneCluster$sysid_main_class_type == p,])[1]
    # get total number of genes that have phenotype p
    m <- dim(totalPhenotypes[totalPhenotypes$sysid_main_class_type==p,])[1]
    # get tota number of genes that dont have phenotype p
    n <- dim(totalPhenotypes[totalPhenotypes$sysid_main_class_type!=p,])[1]
    # get the number of phenotypes known for the cluster mebers
    k <- dim(oneCluster)[1]
    
    pValue <- phyper(q,m,n,k, lower.tail=FALSE)
    #if(pValue < 0.05)
    resultsData[x,p] <- pValue
    
    #print(pValue)    
    if (is.nan(pValue)){
      print(paste('q: ', q))
      print(paste('m: ', m))
      print(paste('n: ', n))
      print(paste('k: ', k))
      nanVector <- c(nanVector,pValue)
      
    }    
  }
}
write.csv(resultsData,  file='/outputData/enchrichment.csv')


clusterOutput <- clustersOfInterest[with(clustersOfInterest, order(cluster_id)),]

write.csv(clusterOutput,  file='/outputData/cluster.csv')





