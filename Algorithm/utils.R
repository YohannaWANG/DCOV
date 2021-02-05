"Author: Yuhao Wang
 Date:   2021-01-28
 Description: 
          - Utils for chain graph struture analysis;
          - And chain graph structure learning baseline algorithms.
          - Implementation for 'Identifiability of chain graph with equal determinant' " 

#library("ggm",lib.loc="F:/Tools/Anaconda3/envs/Bayesian/Lib/R/library")
#library("pcalg",lib.loc="C:/Program Files/R/R-4.0.3/library")
#library("lcd",lib.loc="C:/Program Files/R/R-4.0.3/library")
#Sys.getenv("C:/Program Files/R/R-4.0.3/library")
#Sys.which("stats.dll")
#library("ggm", lib.loc = '/home/yohanna/R/x86_64-pc-linux-gnu-library/4.0')


library(ggm)
library(lcd)
library(pcalg)
source("AMPCGs2019.R")

"TODO:  
1. Use python code generate chain graph adjacency matrix **Adj** and
 simulation data **cg.data**
2. Then return the learnd results **Adj_learned**
3. We can also return the performance evaluation results **results** through comp.cg() function
NOTES: comp.cg() function will return the metrics of [TP, FN, FP, TN, TPR, FPR, SHD]
"


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
# Load "DAG.csv" file. 3000 random samples of the DAG above
cg.data<-read.csv("DAG.csv")
#plot "dag" from the R package lcd
draw(dag)
#check whether "dag" is a chain graph 
is.chaingraph(dag)


"
Function: Return chain components
  Input:  Adj matrix
  Output: Chain components and node number within each chain components
"
chain_comp <- function(cg.data, dag){
  return(is.chaingraph(dag))
}

calculate_mgcv3 <- function(){
  library("mgcv")
  x <- read.csv("DAG.csv")
  ancestors <- unlist(lapply(1, as.numeric))
  current_node <- unlist(lapply(2, as.numeric))
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  mgcvformula = "x[,current_node] ~ "
  for(a in ancestors){
    mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps', sp=0.6) + ")
  }
  mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
  mgcvformula = as.formula(mgcvformula)
  
  b1 = mgcv::gam(mgcvformula, data = x)
  fit.gam = predict(b1)
  
  return(fit.gam)
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


calculate_mgcv2 <- function(){ 
  #source('NPVAR.R')
  #source('utils_NPVAR.R')
  library("mgcv")
  x <- read.csv("DAG.csv")
  ancestors <- unlist(lapply(1, as.numeric))
  current_node <- unlist(lapply(2, as.numeric))
  #ancestors <- lapply(ancestors, as.integer)
  #ancestors <- as.numeric(unlist(ancestors)) + 1
  #current_node <- lapply((current_node + 1), as.integer)
  #current_node  <- as.numeric(unlist(current_node))

  ### Compute nonparametric regression using GAM / mgcv package
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  #colnames(x) <- c(letters[1: ncol(x)])

  mgcvformula = "x[,current_node] ~ "
  for(a in ancestors){
    mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps', sp=0.6) + ")
  }
  mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
  mgcvformula = as.formula(mgcvformula)

  b1 = mgcv::gam(mgcvformula, data = x)

  fit.gam = predict(b1)
  
  return(fit.gam)
}

calculate_mgcv <- function(x, ancestors){ 
  ### Compute nonparametric regression using GAM / mgcv package
  if(!is.data.frame(x)){
    x = data.frame(x)
  }
  print(class(ancestors))

  n = nrow(x)
  p = ncol(x)
  node.index = 1:p
  names(node.index) = names(x)
  
  mgcvformula = "x[,current_node] ~ "
  for(a in ancestors){
    mgcvformula = paste0(mgcvformula, "s(x[,", a, "], bs='ps', sp=0.6) + ")
  }
  mgcvformula = substr(mgcvformula, start = 1, stop = nchar(mgcvformula) - 3)
  mgcvformula = as.formula(mgcvformula)
  
  b1 = mgcv::gam(mgcvformula, data = x)
  fit.gam = predict(b1)
  
  condvar.gam = var(x[,current_node]) - var(fit.gam)
  return(condvar.gam)
}

"NOTES: Algorithms 1-3 are for chain graph structure learning"
"Algoithm implementation code is mainly from the paper:
         2020 - AMP Chain Graphs: Minimal Separators and Structure Learning Algorithms "
"Link: https://github.com/majavid/AMPCGs2019"


"Baseline algorithm 1: LCD-like algorithm"
baseline_lcdlike <- function(cg.data, dag){
  require("ggm")
  require("pcalg")
  require("lcd")
  source("AMPCGs2019.R")
  # Learn the chain graph structure via the LCD-like algorithm
  #colnames(cg.data) <- c("a","b","c","d","e","f")
  colnames(cg.data) <- c(letters[1: ncol(cg.data)])
  row.names(cg.data) <- 1 : nrow(cg.data) 
  colnames(dag) <- c(letters[1: ncol(dag)])
  row.names(dag) <- c(letters[1: ncol(dag)])
  ampcg <- learn.original.amp.normLCD(cg.data, p.value=0.05)
  ampcg <- ampcg[nrow(ampcg):1, ncol(ampcg):1]
  #compare the learned CG to the true CG 
  results <- comp.cgs(dag,ampcg) 
  return(results)
}
"Baseline algorithm 1: LCD-like algorithm"
baseline_lcdlike2 <- function(cg.data, dag){
  require("ggm")
  require("pcalg")
  require("lcd")
  source("AMPCGs2019.R")
  # Learn the chain graph structure via the LCD-like algorithm
  colnames(cg.data) <- c(letters[1: ncol(cg.data)])
  row.names(cg.data) <- 1 : nrow(cg.data) 
  colnames(dag) <- c(letters[1: ncol(dag)])
  row.names(dag) <- c(letters[1: ncol(dag)])
  ampcg <- learn.original.amp.normLCD(cg.data, p.value=0.05)
  #ampcg <- ampcg[nrow(ampcg):1, ncol(ampcg):1]
  #compare the learned CG to the true CG 
  results <- comp.cgs(dag,ampcg) 
  return(ampcg)
}



"Baseline algorithm 2: LDCG algorithm"
baseline_ldcg <- function(cg.data, dag){
  require("ggm")
  require("pcalg")
  require("lcd")
  source("AMPCGs2019.R")
  # Learn the chain graph structure via the LCD-like algorithm
  colnames(cg.data) <- c(letters[1: ncol(cg.data)])
  row.names(cg.data) <- 1 : nrow(cg.data) 
  
  ampcg<-learn.original.amp.normLCD(cg.data,p.value=0.05)
  ampcg<-ampcg[nrow(ampcg):1, ncol(ampcg):1]
  #Learn the largest deflagged graph (LDCG) 
  ldcg<-Largest_DeflaggedAMPCG(ampcg)
  # plot the learned LDCG
  #draw(ldcg)
  # compare the learned LDCG to the true LDCG (dag)
  #results <- comp.cgs(dag,ldcg)  
  return(ldcg)
}

baseline_pclike <-function(cg.data, dag){
  source("AMPCGs2019.R")
  colnames(cg.data) <- c(letters[1: ncol(cg.data)])
  row.names(cg.data) <- 1 : nrow(cg.data) 
  ampcg<-learn.amp.normPC(cg.data, p.value = 0.05, method="stable")
  return(ampcg)
}


" Baseline algorithm 2: LCD algorithm"
baseline_lcd <- function(cg.data, dag){
  require("lcd")
  colnames(cg.data) <- c(letters[1: ncol(cg.data)])
  row.names(cg.data) <- 1 : nrow(cg.data) 
  n <- nrow(cg.data)
  p.value <- .05
  tgug <- naive.getug.norm(cg.data, p.value)
  tg.jtree <- ug.to.jtree(tgug)
  tg.pat <- learn.mec.norm(tg.jtree, cov(cg.data), n, p.value, "CG")
  #df4 <- tg.pat[order(nrow(tg.pat):1) ,order(ncol(tg.pat):1)]
  #comp.skel(skeleton(toy.graph), skeleton(df4))
  #comp.pat(pattern(toy.graph), tg.pat)
  return(tg.pat)
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
  #rownames(amat) <- c("a","b","c","d","e","f")
  #colnames(amat) <- c("a","b","c","d","e","f")
  # Compare the learned AMAT to the true LDCG(dag)
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
  # Reformat row and column name of amat (To enable com.cgs function)
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


