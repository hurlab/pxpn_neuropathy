###########################
#Edge Jaccard 
#Takes a list of graphs
#returns a Jaccard matrix
#of similarity in edge sets
###########################
library(igraph)

jaccard_edge<-function(g1, g2){
  return(length(E(igraph::intersection(g1, g2)))/length(E(igraph::union(g1, g2))))
}

jaccard_edge_matrix <- function(a, b=a){
  matriz = matrix(nrow = length(a), ncol = length(b), dimnames = list(names(a), names(b)))
  for(i in seq_along(a)){
    for(j in seq_along(b)){
      matriz[i,j]<-jaccard_edge(a[[i]], b[[j]])
    }
  }
  return(matriz)
}



