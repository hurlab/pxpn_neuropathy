#jaccard_matrix: calculate Jaccard index between two lists of sets
				#returns a Jaccard matrix 
				#with J between each set in list 


jaccard<-function(a,b)
{
  x<-intersect(a,b)
  y<-union(a,b)
  nx<-length(x)
  ny<-length(y)
  J<-as.numeric(nx/ny)
  print(J)
  return(J)
}

jaccard_matrix<-function(x,y){
  matriz<-matrix(nrow = length(x), ncol = length(y))  
  colnames(matriz)<-names(y)
  rownames(matriz)<-names(x)
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      alfa<-x[[i]]
      beta<-y[[j]]
      matriz[i,j]<-jaccard(alfa,beta)
    }
  }
  return(matriz)  
}

