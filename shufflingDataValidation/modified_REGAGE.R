#############
#REGAGE FUNCTION
#############


library(gage)
library(igraph)

ReGAGE<-function(expmatrix, pathways, cases, controls, qvalue=0.1, setsize = c(10,500)){
  
  ReGAGE_ReSULTS = list()
  
  #Pathway Enrichment
  gage_enrichment = gage(exprs = expmatrix, 
                         gsets = pathways, 
                         ref = controls, 
                         samp = cases,
                         set.size = setsize, 
                         same.dir = FALSE, 
                         compare = "unpaired"
  )
  
  enrichment = as.data.frame(gage_enrichment$greater)
  ReGAGE_ReSULTS$enrichment = enrichment
  
  ##List of EnrichedPathways
  EnrichedPathways = rownames(enrichment)[which(enrichment$q.val<qvalue)]
  ReGAGE_ReSULTS$EnrichedPathways = EnrichedPathways
  
  EnrichedPathways_List = pathways[EnrichedPathways]
  
  #if list of pathways is empty, return an empty graph
  if(length(EnrichedPathways_List)==0){
    ReGAGE_ReSULTS$g = make_empty_graph()
    return(ReGAGE_ReSULTS)
  }
  
  #2) find intersections between EnrichedPathways
  ListIntersectionSets = list()
  for(i in seq_along(EnrichedPathways_List)){
    ListIntersectionSets[[i]] <-lapply(X = tail(EnrichedPathways_List, 
                                                n = length(EnrichedPathways_List) - i), 
                                       FUN = function(x) intersect(x, EnrichedPathways_List[[i]]))
    names(ListIntersectionSets)[i]<-names(EnrichedPathways_List)[i]
  }
  ##List of NONEMPTY IntersectionSets
  ListIntersectionSets_NONEMPTY = Filter(f = length, unlist(ListIntersectionSets, recursive = FALSE))
  
  ##if list of intersection sets is empty, return a graph with all vertices
  
  if(length(ListIntersectionSets_NONEMPTY)==0){
    prueba_g = igraph::add.vertices(graph = make_empty_graph(), 
                                    nv = length(nombres), 
                                    attr = list(name = nombres)
                                    )
    ReGAGE_ReSULTS$g = prueba_g
    return(ReGAGE_ReSULTS)
  }
  
  #3)For every IntersectionSet, hypergeometric (GeneBunch, IntersectionSet)
  gage_enrichment_intersections = gage(exprs = expmatrix, 
                                       gsets = ListIntersectionSets_NONEMPTY, 
                                       ref = controls, 
                                       samp = cases,
                                       set.size = setsize, 
                                       same.dir = FALSE, 
                                       compare = "unpaired"
  )
  
  enrichment_intersections = as.data.frame(gage_enrichment_intersections$greater)
  enrichment_intersections = as.data.frame(enrichment_intersections)
  EnrichedIntersectionSets = rownames(enrichment_intersections)[which(enrichment_intersections$q.val<qvalue)]
  EnrichedIntersectionSets=strsplit(x = EnrichedIntersectionSets, split = "\\.")
  ReGAGE_ReSULTS$enrichment_intersections = enrichment_intersections
  ReGAGE_ReSULTS$EnrichedIntersectionSets = EnrichedIntersectionSets
  
  #if EnrichedIntersectionSets is empty, return a graph with just the enriched nodes
  if(length(EnrichedIntersectionSets)==0){
    prueba_g = igraph::add.vertices(graph = make_empty_graph(), 
                                    nv = length(nombres), 
                                    attr = list(name = nombres)
    )
    ReGAGE_ReSULTS$g = prueba_g
    return(ReGAGE_ReSULTS)
  }
  
  #network generation
  DF <- data.frame(do.call(rbind, EnrichedIntersectionSets)) #contains only significant edges
  DF$weight <- enrichment_intersections[which(enrichment_intersections$q.val<qvalue), "q.val"]
  g = graph_from_edgelist(as.matrix(DF[,1:2]), directed = FALSE)
  E(g)$weight = DF$weight
  
  if(length(setdiff(EnrichedPathways, V(g)$name)) != 0) {
    g = add.vertices(graph = g,
                     nv = length(setdiff(EnrichedPathways, V(g)$name)),
                     name = setdiff(EnrichedPathways, V(g)$name))
    
  }
  
  #add qval as node attribute
  V(g)$qval <-enrichment[match(x = V(g)$name, table = rownames(enrichment)), 4]
  
  ReGAGE_ReSULTS$g = g
  
  return(ReGAGE_ReSULTS)
}
