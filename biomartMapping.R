require('biomaRt')

mappIds <- function(biomart, dataSet, attributes, inputIdType, valuesToMap, uniquerows){
  print('Mapping data via biomart...')
  print(paste('Mart:',biomart)) 
  print(paste('Dataset:',dataSet)) 
  print(paste('Attributes:',attributes)) 
  print(paste('Input Id Type',inputIdType)) 
  
  ensembl <- useMart(biomart,dataset="hsapiens_gene_ensembl")
  return(getBM(attributes, filters = inputIdType, values = valuesToMap, mart = ensembl, uniqueRows=uniquerows))
}