library(gage)
library(tidyverse)
source("modified_REGAGE.R") #use the modified REGAGE version, that returns empty and edgeless networks if there is no enrichment above the q-value threshold

#extract the relevant matrix samples for validation
faketrix = matriz[,c(Cont_SCN_num, dbdb_SCN_num, dbdb_Pio_SCN_num)]


#run 5000 iterations of matrix randomization, 

newNullModels = lapply(X = 1:5000, FUN = function(i){
  print(i)
  #randomize matrix
  nr<-dim(faketrix)[1]
  faketrix2 = faketrix[sample.int(nr),]
  rownames(faketrix2) = rownames(faketrix)
  faketrix2 = as.matrix(faketrix2)
  
  #randomize selected set 
  
  #if "grouping" is to be preserved
  #samplePositionList = list(1:6, 7:12, 13:18)
  #sampledSamples = sample(samplePositionList, 2, replace = FALSE)
  #bb = sampledSamples[[1]]
  #dd = sampledSamples[[2]]
  
  #if groups have any six random samples 
  aa = 1:18
  bb = sample(aa, 6)
  cc = aa[!(aa%in%bb)]
  dd = sample(cc, 6)
  
  #PXPN generation
  fakeGage = ReGAGE(expmatrix = faketrix2, 
                    pathways = reactome_pw, 
                    cases = bb, 
                    controls = dd, 
                    qvalue = 0.05)
  return(fakeGage)
})

newNullResults = lapply(newNullModels, function(l){
  difi = data.frame(no.nodes = length(V(l$g)),
                    no.edges = length(E(l$g)),
                    no.comp  = components(l$g)$no,
                    avgDeg   = mean(degree(l$g)),
                    cc       = transitivity(l$g),
                    avgShPth = average.path.length(l$g)
  )
  return(difi)
})

newDataFrame = plyr::ldply(newNullResults, data.frame)