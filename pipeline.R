####################################################################
#Pathway Crosstalk Perturbation Network Pipeline 
#Diabetic Neuropathy model
####################################################################

#inputs 

## expression matrix 
## grouping information
## pathways

#outputs 

## Pathway Crosstalk Perturbation Networks + communities 

#1) Pathway and crosstalk perturbation > ReGAGEX
#2) Pw Pert NW analysis > NetworkAnalyzer Clone 
#3) Null model generation

#sources and libraries
inputs = scan(file = "inputs.txt", what = "character")
source(file = "sources.R")

# load inputs 

matriz = fread(input = inputs[1], data.table = FALSE)
rownames(matriz) <- matriz[,1]
matriz <- matriz[,-1]
matriz <- matriz[order(rownames(matriz)),]
load(inputs[4])
matriz <- matriz[stringr::str_to_upper(rownames(matriz))%in%features_exp_norm,] #remove non protein coding genes

##groups
groups = fread(inputs[3])
table(groups$Group)

Cont_SCN_num = which(colnames(matriz)%in%groups[Group == "Cont_SCN"]$Sample)
dbdb_SCN_num = which(colnames(matriz)%in%groups[Group == "dbdb_SCN"]$Sample)
dbdb_Pio_SCN_num = which(colnames(matriz)%in%groups[Group == "dbdb-Pio_SCN"]$Sample)

##pathway annotation

reactome_pw = pathways(species = "mmusculus", "reactome")
reactome_pw = convertIdentifiers(x = reactome_pw, to = "symbol")
reactome_pw = lapply(reactome_pw, function(x) nodes(pathwayGraph(x)))


###############################################################################
#1) Pathway and crosstalk perturbation > ReGAGEX
###############################################################################

#ReGAGE, qvalue 0.05



#SCN
SCN_regages = list(alfa = ReGAGE(expmatrix = matriz, 
                                 pathways = reactome_pw, 
                                 cases = dbdb_SCN_num, 
                                 controls = Cont_SCN_num, 
                                 qvalue = 0.05),
                   beta = ReGAGE(expmatrix = matriz, 
                                 pathways = reactome_pw, 
                                 cases = dbdb_Pio_SCN_num, 
                                 controls = dbdb_SCN_num, 
                                 qvalue = 0.05),
                   delta = ReGAGE(expmatrix = matriz, 
                                  pathways = reactome_pw, 
                                  cases = dbdb_Pio_SCN_num, 
                                  controls = Cont_SCN_num, 
                                  qvalue = 0.05)
)

#save(SCN_regages,  file = "...PIO/regages.RData")
#save.image("...PIO/workspace_phase1.RData")


###############################################################################
#2) Pw Pert NW analysis > NetworkAnalyzer Clone 
###############################################################################

SCN_analysis = lapply(SCN_regages, FUN = function(x){
  NetworkAnalyzer(x$g, 
                  directed = FALSE, 
                  skip.betweenness = FALSE, 
                  workaround.betweenness = TRUE)
})

#save(SCN_analysis, file = "...PIO/nw_analysis_regage.RData")

SCN_nws = lapply(SCN_analysis, function(x){
  x$g
})

#save(SCN_nws, file = "...PIO/regage_analyzed_nw_list.RData")
#save.image("...PIO/workspace_phase2.RData")

###############################################################################
#3) Null model generation
###############################################################################

#background probability crosstalk network

background_network_reactome = jaccard_matrix(reactome_pw, reactome_pw)
background_network_reactome = graph_from_adjacency_matrix(background_network_reactome, weighted = TRUE, mode = "undirected")
background_network_reactome = simplify(background_network_reactome, remove.multiple = TRUE, remove.loops = TRUE)
background_network_reactome = NetworkAnalyzer(background_network_reactome)
distribution_plot(background_network_reactome$g)

SCN_bootstraps = lapply(X = SCN_nws, FUN = function(g){
  bootstrap_model_4_function(g = g, background = background_network_reactome$g, n = 5000)
})


#########################################
############DONE!########################
#########################################