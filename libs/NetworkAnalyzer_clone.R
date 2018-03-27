#NetworkAnalyzer clone
#Network parameters calculated by the Cytoscape Network Analyzer
#takes Igraph object
#returns list, including igraph (with centrality measures for vertices), 
#average path length, clustering coefficient, Graph diameter
#component distribution
#and Infomap communities

library(igraph)

NetworkAnalyzer = function(g, directed = FALSE, skip.betweenness = FALSE, workaround.betweenness = FALSE){
  if(directed == FALSE){
    V(g)$degree <- degree(g)
    
    if(workaround.betweenness == TRUE && skip.betweenness == FALSE){
      copy = g
      E(copy)$weight = rep(1, length(E(copy)))
      V(g)$betweenness <- betweenness(copy, directed = FALSE)
    }
    
    if(workaround.betweenness == FALSE && skip.betweenness == FALSE){
    V(g)$betweenness <- betweenness(g, directed = FALSE)
    }
    
    V(g)$closeness <- closeness(g)
    
    average_path_length = average.path.length(g, directed = FALSE)
    ClusteringCoefficient = transitivity(g, type = "global")
    GraphDiameter = diameter(g, directed = FALSE)
    
    component = components(g)
    V(g)$component = component$membership
    
    infomap = infomap.community(g)
    V(g)$infomap = infomap$membership
    
    results = list(g = g, 
                   average_path_length = average_path_length, 
                   ClusteringCoefficient = ClusteringCoefficient, 
                   GraphDiameter = GraphDiameter,
                   component = component,
                   infomap = infomap)
  }
  if(directed == TRUE){
    V(g)$degree <- degree(g, mode = "all")
    V(g)$indegree <- degree(g, mode = "in")
    V(g)$outdegree <- degree(g, mode = "out")
    
    if(workaround.betweenness == TRUE && skip.betweenness == FALSE){
      copy = g
      E(copy)$weight = rep(1, length(E(copy)))
      V(g)$betweenness <- betweenness(copy, directed = TRUE)
    }
    
    if(skip.betweenness == FALSE){
    V(g)$betweenness <- betweenness(g, directed = TRUE)
    }
    
    V(g)$closeness <- closeness(g, mode = "all")
    V(g)$closeness_in <- closeness(g, mode = "in")
    V(g)$closeness_out <- closeness(g, mode = "out")
    
    average_path_length = average.path.length(g, directed = TRUE)
    ClusteringCoefficient = transitivity(g, type = "global")
    GraphDiameter = diameter(g, directed = TRUE)
    
    component = components(g)
    V(g)$component = component$membership
    
    infomap = infomap.community(g)
    V(g)$infomap = infomap$membership
    
    results = list(g = g, 
                   average_path_length = average_path_length, 
                   ClusteringCoefficient = ClusteringCoefficient, 
                   GraphDiameter = GraphDiameter,
                   component = component,
                   infomap = infomap)
  }
  return(results)
}



