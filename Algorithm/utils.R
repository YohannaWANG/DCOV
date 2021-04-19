"Author: 
 Date:   2021-01-28
 Description: 
          - Utils for chain graph struture analysis;
          - And chain graph structure learning baseline algorithms.
          - Implementation for 'Identifiability of chain graph with equal determinant' " 

library(ggm)
# library(lcd)
library(pcalg)
# source("AMPCGs2019.R")

"TODO:  
1. Use python code generate chain graph adjacency matrix **Adj** and
 simulation data **cg.data**
2. Then return the learnd results **Adj_learned**
3. We can also return the performance evaluation results **results** through comp.cg() function
NOTES: comp.cg() function will return the metrics of [TP, FN, FP, TN, TPR, FPR, SHD]
"

test_exp <- function(){ 
  "Load data and graph"
  dag <- matrix(c( 0, 1, 1, 0, 0, 0,
                    0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0),
                 6, 6, byrow = TRUE)
  
  N <- c("a","b","c","d","e","f")
  dimnames(dag) <- list(N, N) 
  cg.data<-read.csv("DAG.csv")
  draw(dag)
  #check whether "dag" is a chain graph 
  is.chaingraph(dag)
}

"Copyright: NPVAR algorithm - generate simulation data"
"Link: https://github.com/MingGao97/NPVAR"
### Simulated data generation
### Given graph type, source variance, number of nodes, sample size, graph degree, x2
### Return a list consists of data matrix and true graph adjacency matrix
data_simu <- function(graph_type, errvar, d, n, s0, x2 = F){
  if (graph_type == 'MC-SIN') {
    G = markov_chain(d)
    X = sampleFromSin(G, n, errvar)
    if (x2) X2 = sampleFromSin(G, n, errvar)
  } else if (graph_type == 'MC-GP') {
    G = markov_chain(d)
    if (x2) {
      X = sampleDataFromG(2 * n, G, errvar = errvar)
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar)
    }
  } else if (graph_type == 'ER-AGP') {
    if ((d==5) & (s0>1)) {
      G = as.matrix(sparsebnUtils::random.graph(d, 9))
    } else {
      G = as.matrix(sparsebnUtils::random.graph(d, s0*d))
    }
    if (x2) {
      X = sampleDataFromG(2 * n, G, errvar = errvar)
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar)
    }
  } else if (graph_type == 'ER-SIN') {
    if ((d==5) & (s0>1)) {
      G = as.matrix(sparsebnUtils::random.graph(d, 9))
    } else {
      G = as.matrix(sparsebnUtils::random.graph(d, s0*d))
    }
    X = sampleFromSin(G, n, errvar)
    if (x2) X2 = sampleFromSin(G, n, errvar)
  } else if (graph_type == 'ER-NGP') {
    if ((d==5) & (s0>1)) {
      G = as.matrix(sparsebnUtils::random.graph(d, 9))
    } else {
      G = as.matrix(sparsebnUtils::random.graph(d, s0*d))
    }
    if (x2) {
      X = sampleDataFromG(2*n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
    }
  } else if (graph_type == 'SF-AGP') {
    G = as_adjacency_matrix(sample_pa(d, m = s0),sparse = F)
    if (x2) {
      X = sampleDataFromG(2 * n, G, errvar = errvar)
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar)
    }
  } else if (graph_type == 'SF-SIN') {
    G = as_adjacency_matrix(sample_pa(d, m = s0),sparse = F)
    X = sampleFromSin(G, n, errvar)
    if (x2) X2 = sampleFromSin(G, n, errvar)
  } else if (graph_type == 'SF-NGP') {
    G = as_adjacency_matrix(sample_pa(d, m = s0),sparse = F)
    if (x2) {
      X = sampleDataFromG(2*n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
      X2 = X[(n+1):(2*n),]
      X = X[1:n,]
    } else {
      X = sampleDataFromG(n, G, errvar = errvar, parsFuncType=list(B=randomB(G),kap=0.01,sigmax=1,sigmay=1,output=FALSE))    
    }
  }
  if (x2) {
    return(list(X=X, G=G, X2=X2))
  } else {
    return(list(X=X, G=G))
  }
}

### Sample from SIN model given adjacency matrix, sample size, source variance
sampleFromSin = function(G, n, errvar = 0.5){
  p = dim(G)[2]
  X = matrix(NA,n,p)
  causOrder = computeCausOrder(G)
  for (node in causOrder) {
    paOfNode = which(G[,node] == 1)
    if(length(paOfNode) == 0){
      X[,node] = rnorm(n, 0, sqrt(errvar))
    }else if(length(paOfNode) == 1){
      X[,node] = sin(X[,paOfNode]) + rnorm(n, 0, sqrt(errvar))
    }else{
      X[,node] = apply(sin(X[,paOfNode]), 1, sum) + rnorm(n, 0, sqrt(errvar))
    }
  }
  return(X)
}

computeCausOrder <- function(G)
  # Copyright (c) 2013  Jonas Peters  [peters@stat.math.ethz.ch]
  # All rights reserved.  See the file COPYING for license terms.
{
  p <- dim(G)[2]
  remaining <- 1:p
  causOrder <- rep(NA,p)
  for(i in 1:(p-1))
  {
    root <- min(which(colSums(G) == 0))
    causOrder[i] <- remaining[root]
    remaining <- remaining[-root]
    G <- G[-root,-root]
  }
  causOrder[p] <- remaining[1]
  return(causOrder)
}


npvar_dag_data <- function(d, n){
  "Generate DAG data based on NPVAR algorithm"
  library("np")
  library("mgcv")
  #source('NPVAR.R')
  #source('utils.R')
  "Generate simulation data"
  data = data_simu(graph_type = 'ER-SIN', errvar = 0.5, d, n, s0 = 1, x2 = T)
  "Extract synthetic data and graph from above"
  X <- data$X
  G <- data$G
  X2 <- data$X2
  ### Naively recover ordering node one by one
  result1 <- NPVAR(X)
  return(list(X, G, result1))
}

"
Function: Return chain components
  Input:  Adj matrix
  Output: Chain components and node number within each chain components
"
chain_comp <- function(cg.data, dag){
  return(is.chaingraph(dag))
}


calculate_np <- function(x, ancestors, current_node){
  library("mgcv")
  #x <- read.csv("DAG.csv")
  ancestors <- unlist(lapply(ancestors, as.numeric)) + 1
  current_node <- unlist(lapply(current_node, as.numeric)) + 1
  
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  
  npformula = "x[,current_node] ~ "
  for(a in ancestors){
    npformula = paste0(npformula, "x[,", a, "] + ")
  }
  npformula = substr(npformula, start = 1, stop = nchar(npformula) - 3)
  npformula = as.formula(npformula)
  
  bw.all = np::npregbw(formula = npformula,
                   regtype = "ll",
                   bwmethod = "cv.aic",
                   data = x)
  model.np = np::npreg(bws = bw.all)
  fit.np = predict(model.np,
                   data = x,
                   newdata = x)
  return(fit.np)
}


"MGCV algorithm:
 Notes: already finished debug"
calculate_mgcv <- function(x, ancestors, current_node){ 
  library("mgcv")
  #x <- read.csv("DAG.csv")
  ancestors <- unlist(lapply(ancestors, as.numeric)) + 1
  currentNode <- unlist(lapply(current_node, as.numeric)) + 1
  ### Compute nonparametric regression using GAM / mgcv package
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  
  mgcvformula = "x[,currentNode] ~ "
  for(a in ancestors){
    as.numeric(unlist(a))
    mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps') + ")
  }
  mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
  
  mgcvformula = as.formula(mgcvformula)
  b1 = mgcv::gam(mgcvformula, data = x,  sp=0.6)
  fit.gam = predict(b1)
  return(fit.gam)
}


"
Function: Prune convert topological order into adjacency matrix
"
### Test if there is violation of estimated ordering under the true adjacency matrix
### Return whether ordering is correct and count of violations
test_order <- function(topo_order, adj){
  edges = which(adj == 1, arr.ind = T)
  order_Index = order(topo_order)
  count = 0
  for (i in 1:nrow(edges)) {
    if (order_Index[edges[i,1]] > order_Index[edges[i,2]]) {
      count = count + 1
    }
  }
  return(list(right = as.numeric(count==0), count = count))
}

### Estimate adjacency matrix using estimated ordering and significance given by GAM
prune <- function(x, est_order, cutoff = 0.001){
  library("mgcv")
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  p = dim(x)[2]
  adj = matrix(0, p, p)
  est_order <- unlist(lapply(est_order, as.numeric)) + 1
  
  for (i in 2:p) {
    node = est_order[i] 
    ancestors = est_order[1:(i-1)]

    mgcvformula = "x[,node] ~ "
    for(a in ancestors){
      mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps') + ")
    }
    mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
    mgcvformula = as.formula(mgcvformula)
    print(mgcvformula)
    mod = mgcv::gam(mgcvformula, data=x, sp=0.6)

    parents = ancestors[summary(mod)$s.pv < cutoff]
    adj[parents, node] = 1
  }
  return(adj)
}

 
"Notes: Algorthms below are DAG structure learning algorithms"

"2020-A polynomial-time algorithm for learning nonparametric causal graphs-NeurIPS 2020"

baseline_NPVAR <- function(cg.data, dag){

  "TODO"
}

"Baseline algorithm 3: GES"
baseline_ges <- function(cg.data, dag){
  library("pcalg")
  score <- new("GaussL0penObsScore", cg.data, lambda = 0.5*log(nrow(cg.data)))
  ## Estimate the essential graph
  ges.fit <- ges(score,phase = c("forward","backward"), iterate = FALSE)
  amat<-wgtMatrix(ges.fit$essgraph) #essgraph (cpdag), repr (An object of a class from ParDAG)
  amat[amat != 0] <- 1
  # Reformat row and column name of amat (To enable com.cgs function)
  return(amat)
}


"Baseline algorithm 4: PC"
baseline_pc <- function(cg.data, dag){
  require("pcalg")
  require("bnlearn")
  require("CompareCausalNetworks")
  methods <- c("pc")
  for (method in methods){ 
    Ahat <- getParents(cg.data, environment=NULL, intervention=NULL, method=method, alpha=0.01, pointConf = TRUE)
  }
  return(Ahat)
}


"Baseline algorithm 5: RFCI"
baseline_rfci <- function(cg.data, dag){
  rfci_cg <- getParents(cg.data, environment=NULL, intervention=NULL, method="rfci", alpha=0.01, pointConf = TRUE)
  # Reformat row and column name of amat (To enable com.cgs function)
  return(rfci_cg)
}

"Baseline algorithm 6: GIES- Greedy Intervention Equaivalence Search"
baseline_gies <- function(cg.data, dag){
  require("pcalg")
  score <- new("GaussL0penObsScore", cg.data, lambda = 0.5*log(nrow(cg.data)))
  gies.fit <- gies(score)
  # plot(gies.fit$essgraph, main = "") ; box(col="gray")
  gies_cg<-wgtMatrix(gies.fit$essgraph) #essgraph (cpdag), repr (An object of a class from ParDAG)
  gies_cg[gies_cg != 0] <- 1
  return(gies_cg)
}

"Baseline algorithm 7: ARGES - Hybrid method for causal discovery
NOTES: This algorithm restrict the search space of GES to subgraphs of a skeleton or 
conditional independence graph (CIG) to estimated in advance. The 'hybrid' method 
means the combination of constrant-based and score-based method."
baseline_arges <- function(cg.data, dag){
  require("pcalg")
  #data <- data.matrix(data, rownames.force = NA)
  ages.fit <- ages(data = cg.data)
  #plot(ages.fit$essgraph, main="Estimated APDAG with AGES")
  arges_cg<-wgtMatrix(ages.fit$essgraph)
  arges_cg[arges_cg != 0] <- 1
  return(amat)  
  
}



