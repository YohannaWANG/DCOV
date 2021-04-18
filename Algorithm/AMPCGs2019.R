#############################################################################################
######################Algorithm 2: ChordlessCycles(G)########################################
#############################################################################################
`ChordlessCycles`<-function(amat){#amat: adjacency matrix of an undirected graph
  counter<-0
  blocked<-c()
  Adj<-list()
  label<-c()
  ccycles<-list()
  ############################################################################################
  ###########Algorithm 3: DegreeLabeling(G)###################################################
  ############################################################################################
  `DegreeLabeling`<-function(amat){
    degree<-c()
    color<-c()
    valency<-c()
    # Adj<-list()
    # label<-c()
    n<-nrow(amat)
    for (i in 1:n) {
      neighbor<-c()
      degree[i]<-valency[i]<-0
      blocked[i]<<-0
      #color[i]=0 means color[i]=white#
      color[i]<-0
      for (j in 1:nrow(amat)) {
        if(amat[i,j]==1){
          degree[i]<-degree[i]+1
          valency[i]<-degree[i]
          neighbor<-c(neighbor,j)
        }
      }
      Adj[[i]]<<-neighbor
    }
    min_degree<-v<-0
    for (i in 1:n) {
      min_degree<-n
      for (j in 1:n) {
        if((color[j]==0) && (degree[j]<min_degree)){
          v<-j
          min_degree<-degree[j]
        }
      }
      label[v]<<-i
      #color[i]=1 means color[i]=black#
      color[v]<-1
      for (k in 1:n) {
        if((amat[v,k]==1) && (color[k]==0)){
          degree[k]<-degree[k]-1
        }
      }
    }
    return(list(label=label,valency=valency,Adj=Adj))
  }
  ##########################################################################################
  ################Algorithm 4: Triplets(G)##################################################
  ##########################################################################################
  `Triplets`<-function(amat){
    triplet<-list()
    n<-nrow(amat)
    valency<-DegreeLabeling(amat)$valency
    Adj<-DegreeLabeling(amat)$Adj
    label<-DegreeLabeling(amat)$label
    t<-1
    c<-1
    for (i in 1:n) {
      #Generate all triplets on form <x, u, y>.
      if(valency[i]>=2){
        sub2Adj<-combn(Adj[[i]],2)
        for (j in 1:ncol(sub2Adj)) {
          if(label[i]<label[sub2Adj[1,j]] && label[sub2Adj[1,j]]<label[sub2Adj[2,j]]){
            if(amat[sub2Adj[1,j],sub2Adj[2,j]]==1){
              ccycles[[c]]<<-c(sub2Adj[1,j],i,sub2Adj[2,j])
              c<-c+1
            }else{
              triplet[[t]]<-c(sub2Adj[1,j],i,sub2Adj[2,j])
              t<-t+1
            }
          }
          if(label[i]<label[sub2Adj[2,j]] && label[sub2Adj[2,j]]<label[sub2Adj[1,j]]){
            if(amat[sub2Adj[1,j],sub2Adj[2,j]]==1){
              ccycles[[c]]<<-c(sub2Adj[2,j],i,sub2Adj[1,j])
              c<-c+1
            }else{
              triplet[[t]]<-c(sub2Adj[2,j],i,sub2Adj[1,j])
              t<-t+1
            }
          }
        }
      } 
    }
    counter<<-length(ccycles)+1
    return(list(triplet=triplet,ccycles=ccycles))
  }
  
  ############################################################################################
  ################Algorithm 6: BlockNeighbors(v, blocked)#####################################
  ############################################################################################
  `BlockNeighbors`<-function(v){
    for (i in 1:length(Adj[[v]])) {
      blocked[Adj[[v]][i]]<<-blocked[Adj[[v]][i]]+1
    }
  }
  
  ############################################################################################
  ########################Algorithm 7: UnblockNeighbors(v, blocked)###########################
  ############################################################################################
  `UnblockNeighbors`<-function(v){
    for (i in 1:length(Adj[[v]])){
      if(blocked[Adj[[v]][i]]>0){
        blocked[Adj[[v]][i]]<<-blocked[Adj[[v]][i]]-1
      }
    }
  }
  
  ############################################################################################
  #######################Algorithm 5: CC_Visit(p, C, key, blocked)############################
  ############################################################################################
  `CC_Visit`<-function(cpath,ccycles){
    key<-label[cpath[2]]
    BlockNeighbors(cpath[length(cpath)])
    for (i in Adj[[cpath[length(cpath)]]]) {
      if(label[i]>key && blocked[i]==1){
        cpath1<-c(cpath,i)
        if(cpath[1] %in% Adj[[i]]){
          ccycles[[counter]]<<-cpath1
          counter<<-counter+1
        }else{
          CC_Visit(cpath1,ccycles)
        }
      }
    }
    UnblockNeighbors(cpath[length(cpath)])
  }
  ###################################################Algorithm 2: continue...#################
  G<-DegreeLabeling(amat)
  triplet<-Triplets(amat)$triplet
  t<-1
  while (length(triplet)>0) {
    cpath<-triplet[[t]]
    BlockNeighbors(cpath[2])
    ccycles[counter]<-CC_Visit(cpath,ccycles)
    UnblockNeighbors(cpath[2])
    triplet[[t]]<-NULL
  }
  return(ccycles)
}
################################################################################################################
`comp.cgs` <- function(source, target)
{ ###Compares the (learned) largest deflaged chain graph to the (supposed) true largest deflaged chain graph. The two patterns should
  ###have the same vertex set in order for the function to return a meaningful result.
  vset <- colnames(source)
  
  ###true directed edges
  truearr <- which(source - t(source) == 1)
  
  ###true undirected edges
  trueline<-intersect(which(source == 1), which(source == t(source)))
  
  target <- target[vset, vset]
  `skelet` <- function(amat)
  {
    0 + (amat + t(amat) > 0)
  }
  
  ###true learned skelet
  trueskel<-which(skelet(source)+skelet(target)==2)
  
  truepos<-length(trueskel)/2
  ###sensitivity, recall, hit rate, or true positive rate (TPR): TPR=TP/P=TP/TP+FN
  ####condition positive (P): the number of real positive cases
  P<-length(truearr)+(length(trueline)/2)
  TPR<-truepos/P
  total<-choose(length(vset),2)
  
  falsepos<-(length(which(skelet(target)!=0))/2)-truepos
  ####fall-out or false positive rate (FPR): FPR=FP/N=FP/FP+TN
  ####condition negative (N): the number of real negative cases
  N<-total-P
  FPR<-falsepos/N
  
  ####accuracy (ACC): ACC=TP+TN/P+N
  ####N=TN+FP
  trueneg<-N-falsepos
  ACC<-(truepos+trueneg)/(P+N)
  
  ####False Negative: P=TP+FN
  falseneg<-P-truepos
  
  ## computing structural Hamming distance
  #Structural Hamming distance is defined as the total number of operations needed to convert one
  # graph to the other. Each operation must be one of the following: (1) add or delete an
  # edge, or (2) add, remove or reverse an orientation of an edge.
  r=nrow(source)
  shd<-0
  for (i in 1:(r-1)) {
    for (j in (i+1):r) {
      if((source[i,j]!=0 || source[j,i]!=0)&&(target[i,j]==0)&&(target[j,i]==0)){
        ###missing edge in the target: add appropriate edge
        shd<-shd+1
      }
      if((target[i,j]!=0 || target[j,i]!=0)&&(source[i,j]==0)&&(source[j,i]==0)){
        ###extra edge in the target: remove it
        shd<-shd+1
      }
      if((source[i,j]+source[j,i]==1)&&(target[i,j]+target[j,i]==2)){
        ###there is a directed edge in the source graph, but the corresponding edge in target is undirected: add proper orientation 
        shd<-shd+1
      }
      if((source[i,j]+source[j,i]==1)&&(target[i,j]+target[j,i]==1)&&(source[i,j]!=target[i,j])){
        ###-->/<-- in the source graph, but <--/--> in the target, respectively: reverse the orientation 
        shd<-shd+1
      }
      if((source[i,j]+source[j,i]==2) && (target[i,j]+target[j,i]==1)){
        ###there is an undirected edge in the source, but the corresponding edge in target is directed -->/<--: remove the orientation
        shd<-shd+1
      }
    }
  }
  return(list(TP = truepos,
              FN = falseneg,
              FP = falsepos,
              TN = trueneg,
              TPR = TPR,
              FPR = FPR,
              ACC = ACC,
              SHD = shd))
}
##################################################################################################################################
Gskeleton <- function(suffStat, indepTest, alpha, labels, p,
                      method = c("stable", "original", "stable.fast"), m.max = Inf,
                      fixedGaps = NULL, fixedEdges = NULL,
                      NAdelete = TRUE, numCores = 1, verbose = FALSE)
{
  ## Purpose: Perform undirected part of PC-Algorithm, i.e.,
  ## estimate skeleton of AMPCG given data
  ## Order-independent version! NEU
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat: List containing all necessary elements for the conditional
  ##             independence decisions in the function "indepTest".
  ## - indepTest: predefined function for testing conditional independence
  ## - alpha: Significance level of individual partial correlation tests
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - numCores: number of cores to be used for calculation if 
  ##   method = "stable.fast"
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G, sepset, pMax, ord, n.edgetests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 09.12.2009
  ## Modification: Diego Colombo; Martin Maechler; Alain Hauser
  
  ## x,y,S konstruieren
  ##-   tst <- try(indepTest(x,y,S, obj))
  ##-   if(inherits(tst, "try-error"))
  ##-     stop("the 'indepTest' function does not work correctly with 'obj'")
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p))
      p <- length(labels)
    else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    ## Don't want message, in case this is called e.g. from fciPlus():
    ## else
    ##   ###message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  ## C++ version still has problems under Windows; will have to check why
  ##  if (method == "stable.fast" && .Platform$OS.type == "windows") {
  ##    method <- "stable"
  ##    warning("Method 'stable.fast' is not available under Windows; using 'stable' instead.")
  ##  }
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
    ###message("G=")
    ###print(G)
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  
  diag(G) <- FALSE
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)) )
    stop("fixedEdges must be symmetric")
  
  ## Check number of cores
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && numCores > 0)
  if (numCores > 1 && method != "stable.fast") {
    warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
  }
  if (method == "stable.fast") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose),
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
      NAdelete = NAdelete,
      numCores = numCores)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options);
    G <- res$amat
    ## sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1L
  }
  else {
    ## Original R version
    
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)# final length = max { ord}
    ###message("class(n.edgetests)= ",class(n.edgetests))
    ###print(n.edgetests)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      #print(ind)
      remEdges <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,"\n",sep = "")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G.l <- split(G, gl(p,p))
      }
      for (i in 1:remEdges) {
        if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if(method == "stable") G.l[[x]] else G[,x]
          ###message("nbrsBool ",V[ind[i, 1]], "= x, y= ",V[ind[i, 2]])
          ###print(nbrsBool)
          nbrsBool[y] <- FALSE
          ###find adjacent variables
          nbrs <- seq_p[nbrsBool]
          copyOfnbrs<-nbrs
          ###message("nbrs= ")
          ###print(nbrs)
          LNbs<-length(nbrs)
          ###also find adjacents of adjacents
          if(LNbs>0){
            copyOfnbrs<-nbrs
            for (k in 1:LNbs) {
              tempNbsBool<-if(method == "stable") G.l[[nbrs[k]]] else G[,nbrs[k]]
              tempNbsBool[x]<-tempNbsBool[y]<-FALSE
              copyOfnbrs<-union(copyOfnbrs,seq_p[tempNbsBool])
            }
            #nbrs<-sort(copyOfnbrs)
            nbrs<-copyOfnbrs
          }
          ###message("nbrs= ")
          ###print(nbrs)
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            ###message("S=")
            ###print(S)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if(pval >= alpha) { # independent
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            } ## repeat
          }
        }
      }# for( i )
      ord <- ord + 1L
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }
  
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }
  
  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}## end{ Gskeleton }
####################################################################################################################
`learn.amp.multinomPC`<-function(data,p.value,method ="stable",...){
  #### First, load the package pcalg and the data set. ####
  V <- colnames(data)
  #covariance<-cov(data)
  #print(covariance)
  n<-nrow(data)
  df <- data.frame(lapply(data, function(x) if(is.logical(x)) { 
    return(as.factor(x))
  } else {  
    return(x) 
  }
  ), stringsAsFactors=FALSE)
  nlevels<-sapply(df, nlevels)
  dm<-data.matrix(df, rownames.force = NA)
  if(any(dm==0)){
    dm<-dm
  }else{
    dm<- dm-1
  }
  #print(dm[1:50,])
  suffStat <- list(dm =dm , nlev = nlevels, adaptDF = FALSE)
  mm <- switch(method,
               stable = 1,
               original = 2,
               stable.fast = 3,
               0)
  if (mm == 0) stop("Invalid method!")
  if(mm==1){
    skel<-Gskeleton(suffStat,
                    indepTest = disCItest, ## (partial correlations)
                    labels =V,alpha=p.value, verbose = FALSE,...)
  }
  if(mm==2){
    skel<-Gskeleton(suffStat,labels =V,method="original",
                    indepTest = disCItest, ## (partial correlations)
                    alpha=p.value, verbose = FALSE,...)
  }
  if(mm==3){
    skel<-Gskeleton(suffStat,labels =V,method="stable.fast",
                    indepTest = disCItest, ## (partial correlations)
                    alpha=p.value, verbose = FALSE,...)
  }
  wmat<-as(skel@graph,"matrix")####gives adjacency matrix of pcalg, which is corresponding to the learned skeleton
  #print(skel@sepset)
  
  #################################################################
  #################################################################
  #finding chordless cycles 
  ccycles<-ChordlessCycles(wmat)
  lcc<-length(ccycles)
  ########################################################################
  ########## Finding possible parent-neighbors nodes for a given node in 
  ########## a chordless cycle ###########################################
  ########################################################################
  `findPNs`<-function(v,ccycle){
    vParentsNeighbors<-c()
    for (i in 1:length(ccycle)) {
      if(wmat[v,ccycle[i]]==5){
        vParentsNeighbors<-c(vParentsNeighbors,ccycle[i])
      }
    }
    return(vParentsNeighbors)
  }
  ##########################################################################################################
  ##########################################################################################################
  ###A function that determines the type of chordless cycles: with or without blocks########################
  ### if the cycle is chordless and without blocks, the function returns 1 o.w., returns 0 #################
  ##########################################################################################################
  ##########################################################################################################
  ccycles_type<-function(vect,amat){
    cu<-cd<-co<-0
    l<-length(vect)
    for (i in 1:(l-1)) {
      if(amat[vect[i],vect[i+1]]==1 &&  amat[vect[i+1],vect[i]]==1){
        cu<-cu+1
      }
    }
    if(amat[vect[l],vect[1]]==1 &&  amat[vect[1],vect[l]]==1){
      cu<-cu+1
    }
    if(cu==l){
      #undirected cycle
      return(1)
    }else{
      return(0)
    }
  }
  ########################################################################################
  ########################################################################################
  ########## Then apply the rules: 1-4 ###################################################
  ########################################################################################
  change<-TRUE
  while(change){
    change<-FALSE
    ########################################################################
    ###### Rule 1 ###### A o--o B o--o C & B\not\in S_AC -> A +-o B o-+ C
    ########################################################################
    for (i in 1:length(V)) {
      #message(V[i],": *")
      bAdjs<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          bAdjs<-c(bAdjs,j)
        }
      }
      ####Then check A and C are not adjacent and B\not\in S_AC
      if(length(bAdjs)>1){
        for (k in 1:(length(bAdjs)-1)) {
          for(l in (k+1):length(bAdjs)){
            aIndex<-bAdjs[k]
            cIndex<-bAdjs[l]
            # message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
            # print(skel@sepset[[aIndex]][[cIndex]])
            # print(skel@sepset[[cIndex]][[aIndex]])
            ##if they are not adjacents and B\not\in S_AC
            if((wmat[i,aIndex]!=5 || wmat[i,cIndex]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && !((i %in% skel@sepset[[aIndex]][[cIndex]]) ||
                                                                                                                 (i %in% skel@sepset[[cIndex]][[aIndex]]))){
              ###then block from B to A and C
              wmat[i,aIndex]<-5
              wmat[i,cIndex]<-5
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 1
    ########################################################################
    ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
    ########################################################################
    for (i in 1:length(V)) {
      #message("change= ",change)
      bAdjs<-c()#parent, neighbor, or children of B
      bParentsOrNeighbors<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          bAdjs<-c(bAdjs,j)
        }
        if(wmat[i,j]==5){
          bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
        }
      }
      #message(i," bParNeighbors: ",bParentsOrNeighbors)
      ####Then check A and C are not adjacent and B\in S_AC
      if(length(bAdjs)>0 && length(bParentsOrNeighbors)>0){
        for (k in 1:length(bParentsOrNeighbors)) {
          for(l in 1:length(bAdjs)){
            aIndex<-bParentsOrNeighbors[k]
            cIndex<-bAdjs[l]
            ##if A,C are not adjacents and B\in S_AC
            if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && ((i %in% skel@sepset[[aIndex]][[cIndex]]) ||
                                                                                                               (i %in% skel@sepset[[cIndex]][[aIndex]]))){
              # message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
              # print(skel@sepset[[aIndex]][[cIndex]])
              # print(skel@sepset[[cIndex]][[aIndex]])
              ###then block from C to B
              wmat[cIndex,i]<-5
              #message("b: ",V[i]," ","neighbor: ",V[nbIndex])
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 2
    ########################################################################
    ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
    ########################################################################
    if(lcc>0){
      for (i in 1:lcc) {
        cc<-ccycles[[i]]
        #print(cc)
        #cc[j] plays the role of A
        for (j in 1:length(cc)) {
          aAdjs<-c()
          for (k in 1:length(cc)) {
            #aAdjs<-c()
            if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
              aAdjs<-c(aAdjs,cc[k])
            }
          }
          #message(cc[j]," as A: ",aAdjs)
          if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
            #aAdjs[1] is a candidate for being B
            bPN<-findPNs(aAdjs[1],cc)
            if(length(bPN)==1){
              wmat[aAdjs[1],cc[j]]<-5
              change<-TRUE
              break
            }
          }
          if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
            #aAdjs[2] is a candidate for being B
            bPN<-findPNs(aAdjs[2],cc)
            if(length(bPN)==1){
              wmat[aAdjs[2],cc[j]]<-5
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 3
    
    ########################################################################
    ###### Rule 4 ##########################################################
    ########################################################################
    for (i in 1:length(V)) {
      aAdjs<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          aAdjs<-c(aAdjs,j)
        }
      }
      if(length(aAdjs)>=3){
        for (k in 1:length(aAdjs)) {
          for (l in 1:length(aAdjs)) {
            cIndex<-aAdjs[k]
            dIndex<-aAdjs[l]
            ##if C,D are not adjacents and A\in S_CD
            if((cIndex!=dIndex) && ((i %in% skel@sepset[[cIndex]][[dIndex]]) ||
                                    (i %in% skel@sepset[[dIndex]][[cIndex]]))){
              for (r in 1:length(aAdjs)) {
                bIndex<-aAdjs[r]
                #and if B to C,D is blocked and B to A is not blocked
                if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
                  #then block B to A i.e., A +- B
                  wmat[bIndex,i]<-5
                  change<-TRUE
                  break
                }
              }
            }
          }
        }
      }
    }##end of rule 4
  }
  # ##################################################################################################
  #################################################################################################### 
  # ################### This part of code can be used to obtain the Essential graph ##################
  # ################### of the Markov equivalence class of G and not just a CG in class ##############
  # ##### Learning marginal AMP chain graphs under faithfulness revisited (Pena & Gomez-Olmedo, 2016)#
  # ##################################################################################################
  # #  Replace every edge A - B in every cycle that is of length greater than three, chordless,####### 
  # #  and without blocks with A +-+ B ###############################################################
  # ##################################################################################################
  #   if(lcc>0){
  #     for (i in 1:lcc) {
  #       c<-ccycles[i]
  #       l<-length(c)
  #       if((length(c)>3) && (ccycles_type(c,wmat)==1)){
  #         for (j in 1:(l-1)) {
  #           wmat[c[j],c[j+1]]<-wmat[c[j+1],c[j]]<-5
  #         }
  #         wmat[c[1],c[l]]<-wmat[c[l],c[1]]<-5
  #       }
  #     }
  #   }
  # ##################################################################################################
  # ##################################################################################################
  # ##  Apply the rules R2-R4 while possible #########################################################
  # ##################################################################################################
  # ##################################################################################################
  #   change<-TRUE
  #   while(change){
  #     change<-FALSE
  #     ########################################################################
  #     ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       #message("change= ",change)
  #       bAdjs<-c()#parent, neighbor, or children of B
  #       bParentsOrNeighbors<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           bAdjs<-c(bAdjs,j)
  #         }
  #         if(wmat[i,j]==5){
  #           bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
  #         }
  #       }
  #       #message(i," bParNeighbors: ",bParentsOrNeighbors)
  #       ####Then check A and C are not adjacent and B\in S_AC
  #       if(length(bAdjs)>0 && length(bParentsOrNeighbors)>0){
  #         for (k in 1:length(bParentsOrNeighbors)) {
  #           for(l in 1:length(bAdjs)){
  #             aIndex<-bParentsOrNeighbors[k]
  #             cIndex<-bAdjs[l]
  #             ##if A,C are not adjacents and B\in S_AC
  #             if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && ((i %in% skel@sepset[[aIndex]][[cIndex]]) ||
  #                                                                                                                (i %in% skel@sepset[[cIndex]][[aIndex]]))){
  #               # message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
  #               # print(skel@sepset[[aIndex]][[cIndex]])
  #               # print(skel@sepset[[cIndex]][[aIndex]])
  #               ###then block from C to B
  #               wmat[cIndex,i]<-5
  #               #message("b: ",V[i]," ","neighbor: ",V[nbIndex])
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 2
  #     ########################################################################
  #     ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
  #     ########################################################################
  #     if(lcc>0){
  #       for (i in 1:lcc) {
  #         cc<-ccycles[[i]]
  #         #print(cc)
  #         #cc[j] plays the role of A
  #         for (j in 1:length(cc)) {
  #           aAdjs<-c()
  #           for (k in 1:length(cc)) {
  #             #aAdjs<-c()
  #             if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
  #               aAdjs<-c(aAdjs,cc[k])
  #             }
  #           }
  #           #message(cc[j]," as A: ",aAdjs)
  #           if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
  #             #aAdjs[1] is a candidate for being B
  #             bPN<-findPNs(aAdjs[1],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[1],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #           if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
  #             #aAdjs[2] is a candidate for being B
  #             bPN<-findPNs(aAdjs[2],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[2],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 3
  #     
  #     ########################################################################
  #     ###### Rule 4 ##########################################################
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       aAdjs<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           aAdjs<-c(aAdjs,j)
  #         }
  #       }
  #       if(length(aAdjs)>=3){
  #         for (k in 1:length(aAdjs)) {
  #           for (l in 1:length(aAdjs)) {
  #             cIndex<-aAdjs[k]
  #             dIndex<-aAdjs[l]
  #             ##if C,D are not adjacents and A\in S_CD
  #             if((cIndex!=dIndex) && ((i %in% skel@sepset[[cIndex]][[dIndex]]) ||
  #                                     (i %in% skel@sepset[[dIndex]][[cIndex]]))){
  #               for (r in 1:length(aAdjs)) {
  #                 bIndex<-aAdjs[r]
  #                 #and if B to C,D is blocked and B to A is not blocked
  #                 if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
  #                   #then block B to A i.e., A +- B
  #                   wmat[bIndex,i]<-5
  #                   change<-TRUE
  #                   break
  #                 }
  #               }
  #             }
  #           }
  #         }
  #       }
  #     }##end of rule 4
  #   } 
  ##################################################################################################
  ################## Replace every edge A +- B with A -> B #########################################
  ################## Replace every edge A +-+ B with A - B #########################################
  ##################################################################################################
  for (i in 1:(length(V)-1)) {
    for (j in (i+1):length(V)) {
      if(wmat[i,j]==5 && wmat[j,i]==5){
        wmat[i,j]<-wmat[j,i]<-1
      }
      if(wmat[i,j]==5 && wmat[j,i]==1){
        wmat[i,j]<-0
      }
      if(wmat[j,i]==5 && wmat[i,j]==1){
        wmat[j,i]<-0
      }
    }
  }
  return(wmat)
}
####################################################################################################################
`learn.amp.multinomLCD`<-function(tgdata,p.value){
  #### First, load the package pcalg and the data set. ####
  V <- colnames(tgdata)
  freq.table<-as.freq.tb(data.matrix(tgdata, rownames.force = NA))
  tgug<-naive.getug.multinom(freq.table,p.value = p.value,method = "mkb")
  tg.jtree <- ug.to.jtree(tgug)
  skel<-learn.skeleton.multinomAMPCGs(tg.jtree, freq.table, p.value, drop = TRUE)
  wmat <- skel$amat#gives adjacency matrix of the learned skeleton via LCD algorithm
  vset <- rownames(wmat)
  sep.pairs <- skel$sep.pairs
  cliques <- tg.jtree@cliques
  n.clique <- length(cliques)
  ###print(wmat)
  ###print(tg.jtree)
  #################################################################
  #################################################################
  #finding chordless cycles 
  ccycles<-ChordlessCycles(wmat)
  lcc<-length(ccycles)
  #finding chordless cycles of length 3
  ccyclesOf3<-list()
  vtemp<-c()
  if(lcc>0){
    for (i in 1:lcc){
      if(length(ccycles[[i]])==3){
        vtemp<-c(vtemp,i)
      }
    }
  }
  if(length(vtemp)>0){
    ccyclesOf3<-ccycles[vtemp]
  }
  lccOf3<-length(ccyclesOf3)
  ########################################################################
  ########## Finding possible parent-neighbors nodes for a given node in 
  ########## a chordless cycle ###########################################
  ########################################################################
  `findPNs`<-function(v,ccycle){
    vParentsNeighbors<-c()
    for (i in 1:length(ccycle)) {
      if(wmat[v,ccycle[i]]==5){
        vParentsNeighbors<-c(vParentsNeighbors,ccycle[i])
      }
    }
    return(vParentsNeighbors)
  }
  # ##########################################################################################################
  # ##########################################################################################################
  # ###A function that determines the type of chordless cycles: with or without blocks########################
  # ### if the cycle is chordless and without blocks, the function returns 1 o.w., returns 0 #################
  # ##########################################################################################################
  # ##########################################################################################################
  # ccycles_type<-function(vect,amat){
  #   cu<-cd<-co<-0
  #   l<-length(vect)
  #   for (i in 1:(l-1)) {
  #     if(amat[vect[i],vect[i+1]]==1 &&  amat[vect[i+1],vect[i]]==1){
  #       cu<-cu+1
  #     }
  #   }
  #   if(amat[vect[l],vect[1]]==1 &&  amat[vect[1],vect[l]]==1){
  #     cu<-cu+1
  #   }
  #   if(cu==l){
  #     #undirected cycle
  #     return(1)
  #   }else{
  #     return(0)
  #   }
  # }
  # ########################################################################################
  # ########## A function that returns the separation set of two given variables ###########
  # ########################################################################################
  `AsepCD`<-function(u,v,w){
    pair <- new("sep.pair", u=vset[u], v=vset[v])
    a <- lapply(sep.pairs, function(x) all.equal(x,pair))
    ###message(a," for ",u," and ", v)
    if(length((which(a==TRUE)))>0){
      sep<-sep.pairs[[which(a==TRUE)[1]]]@s
      if(is.element(vset[w], sep)){
        return(1)
      }else{
        return(0)
      }
    }else{
      return(2)
    }
  }
  
  ########################################################################################
  ########################################################################################
  ########## Then apply the rules: 1-4 ###################################################
  ########################################################################################
  change<-TRUE
  while(change){
    change<-FALSE
    ########################################################################
    ###### Rule 1 ###### A o--o B o--o C & B\not\in S_AC -> A +-o B o-+ C
    ########################################################################
    for(i in 1:n.clique){
      cvset <- cliques[[i]]$vset
      ##print(cvset)
      p <- length(cvset)
      if(p > 2){
        for(j in 1:(p-1)){
          for(k in (j+1):p){
            for(l in 1:p){
              u <- cvset[j]
              v <- cvset[k]
              w <- cvset[l]
              ####Then check A and C are not adjacent and B\not\in S_AC
              if((wmat[u,v]+wmat[v,u] == 0) &&
                 (wmat[w,u] != 5 || wmat[w,v] != 5) &&
                 (wmat[u,w]+wmat[w,u] > 0) && (wmat[w,v]+wmat[v,w] > 0)){
                pair <- new("sep.pair", u=u, v=v)
                a <- lapply(sep.pairs, function(x) all.equal(x,pair))
                ###message(a," for ",u," and ", v)
                sep <- sep.pairs[[which(a==TRUE)[1]]]@s
                ###message(sep," for ",u," and ", v)
                if(!is.element(w, sep)){
                  ###then block from B to A and C
                  wmat[w,u] <-5
                  wmat[w,v] <-5
                  change<-TRUE
                  ##message(w ," is a sudo-collider for ", u ," and ", v," rule 1")
                  ##print(wmat)
                  break
                }
              }#second if
              #count<- count + 1
            }#fourth for loop
          }#third for loop
        }#second for loop
      }#first if
      ###message("counter= ",count)
    }#first for loop
    #end of rule 1
    ########################################################################
    ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
    ########################################################################
    for (i in 1:length(V)) {
      ###message("change= ",change)
      bAdjs<-c()#parent, neighbor, or children of B
      bParentsOrNeighbors<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          bAdjs<-c(bAdjs,j)
        }
        if(wmat[i,j]==5){
          bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
        }
      }
      ###message(i," bParNeighbors: ",bParentsOrNeighbors)
      ####Then check A and C are not adjacent and B\in S_AC
      if(length(bAdjs)>1 && length(bParentsOrNeighbors)>0){
        for (k in 1:length(bParentsOrNeighbors)) {
          for(l in 1:length(bAdjs)){
            aIndex<-bParentsOrNeighbors[k]
            cIndex<-bAdjs[l]
            ##if A,C are not adjacents and B\in S_AC
            if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && (wmat[aIndex,cIndex]+ wmat[cIndex,aIndex]==0) && (AsepCD(aIndex,cIndex,i)==1 ||
                                                                                                            AsepCD(aIndex,cIndex,i)==2)){
              # ##message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
              # ##print(skel@sepset[[aIndex]][[cIndex]])
              # ##print(skel@sepset[[cIndex]][[aIndex]])
              ###then block from C to B
              wmat[cIndex,i]<-5
              ##message("b: ",vset[i]," "," c: ",vset[cIndex]," a: ",vset[aIndex]," rule 2")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 2
    ########################################################################
    ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
    ########################################################################
    if(lcc>0){
      for (i in 1:lcc) {
        cc<-ccycles[[i]]
        ###print(cc)
        #cc[j] plays the role of A
        for (j in 1:length(cc)) {
          aAdjs<-c()
          for (k in 1:length(cc)) {
            #aAdjs<-c()
            if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
              aAdjs<-c(aAdjs,cc[k])
            }
          }
          ###message(cc[j]," as A: ",aAdjs)
          if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
            #aAdjs[1] is a candidate for being B
            bPN<-findPNs(aAdjs[1],cc)
            if(length(bPN)==1){
              wmat[aAdjs[1],cc[j]]<-5
              ##message("aAdjs[1] ",aAdjs[1]," cc[j] ",cc[j], "rule 3a")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
          if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
            #aAdjs[2] is a candidate for being B
            bPN<-findPNs(aAdjs[2],cc)
            if(length(bPN)==1){
              wmat[aAdjs[2],cc[j]]<-5
              ##message("aAdjs[1] ",aAdjs[1]," cc[j] ",cc[j], "rule 3b")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 3
    ###Now, since rule 1, 2, and 3 are faster than rule 4, apply these until no change, then go on to rule 4
    if (change){
      next
    }
    
    ########################################################################
    ###### Rule 4 ##########################################################
    ########################################################################
    
    if(lccOf3>1){
      for (i in 1:((lccOf3)-1)) {
        for (j in (i+1):lccOf3) {
          abIndex<-intersect(ccyclesOf3[[i]],ccyclesOf3[[j]])
          cIndex<-setdiff(ccyclesOf3[[i]],abIndex)
          dIndex<-setdiff(ccyclesOf3[[j]],abIndex)
          if((length(abIndex)==2) && (wmat[cIndex,dIndex]+wmat[dIndex,cIndex]==0)){
            for (turn in 1:2) {
              aIndex <- if(turn == 1) abIndex[1] else abIndex[2]
              bIndex <- if(turn == 1) abIndex[2] else abIndex[1]
              ##if C,D are not adjacents and A\in S_CD and if B to C,D is blocked and B to A is not blocked
              if(wmat[bIndex,aIndex]!=5 && wmat[aIndex,bIndex]!=0 && wmat[bIndex,cIndex]==5 && wmat[bIndex,dIndex]==5 &&
                 (AsepCD(cIndex,dIndex,aIndex)==1 || AsepCD(cIndex,dIndex,aIndex)==2)){
                #then block B to A i.e., A +- B
                wmat[bIndex,aIndex]<-5
                ##message("bIndex ",vset[bIndex]," aIndex ",vset[aIndex]," cIndex ",vset[cIndex]," dIndex ",vset[dIndex])
                ##print(wmat)
                change<-TRUE
                break
              }
            }
          }
        }
      }
    }##end of rule 4
    
    # for (i in 1:length(V)) {
    #   aAdjs<-c()
    #   for (j in 1:length(V)) {
    #     if(wmat[j,i]+wmat[i,j]>0){
    #       aAdjs<-c(aAdjs,j)
    #     }
    #   }
    #   if(length(aAdjs)>=3){
    #     for (k in 1:length(aAdjs)) {
    #       for (l in 1:length(aAdjs)) {
    #         cIndex<-aAdjs[k]
    #         dIndex<-aAdjs[l]
    #         ##if C,D are not adjacents and A\in S_CD
    #         if((cIndex!=dIndex) && !((wmat[i,cIndex]+wmat[cIndex,i]==2 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,dIndex]+wmat[dIndex,i]==2 && wmat[i,cIndex]==5 && wmat[cIndex,i]==1) ||
    #                                 (wmat[i,cIndex]==5 && wmat[cIndex,i]==1 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,cIndex]+wmat[cIndex,i]==10 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,dIndex]+wmat[dIndex,i]==10 && wmat[i,cIndex]==5 && wmat[cIndex,i]==1))){
    #           for (r in 1:length(aAdjs)) {
    #             bIndex<-aAdjs[r]
    #             #and if B to C,D is blocked and B to A is not blocked
    #             if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
    #               #then block B to A i.e., A +- B
    #               wmat[bIndex,i]<-5
    #               ##message("bIndex ",vset[bIndex]," A ",vset[i]," cIndex ",vset[cIndex]," dIndex ",vset[dIndex])
    #               ##print(aAdjs)
    #               ##print(wmat)
    #               change<-TRUE
    #               #break
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }##end of rule 4
  }
  # ##################################################################################################
  #################################################################################################### 
  # ################### This part of code can be used to obtain the Essential graph ##################
  # ################### of the Markov equivalence class of G and not just a CG in class ##############
  # ##### Learning marginal AMP chain graphs under faithfulness revisited (Pena & Gomez-Olmedo, 2016)#
  # ##################################################################################################
  # ##################################################################################################
  # ##################################################################################################
  # #  Replace every edge A - B in every cycle that is of length greater than three, chordless,####### 
  # #  and without blocks with A +-+ B ###############################################################
  # ##################################################################################################
  #   if(lcc>0){
  #     for (i in 1:lcc) {
  #       c<-ccycles[i]
  #       l<-length(c)
  #       if((length(c)>3) && (ccycles_type(c,wmat)==1)){
  #         for (j in 1:(l-1)) {
  #           wmat[c[j],c[j+1]]<-wmat[c[j+1],c[j]]<-5
  #         }
  #         wmat[c[1],c[l]]<-wmat[c[l],c[1]]<-5
  #       }
  #     }
  #   }
  # ##################################################################################################
  # ##################################################################################################
  # ##  Apply the rules R2-R4 while possible #########################################################
  # ##################################################################################################
  # ##################################################################################################
  #   change<-TRUE
  #   while(change){
  #     change<-FALSE
  #     ########################################################################
  #     ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       ###message("change= ",change)
  #       bAdjs<-c()#parent, neighbor, or children of B
  #       bParentsOrNeighbors<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           bAdjs<-c(bAdjs,j)
  #         }
  #         if(wmat[i,j]==5){
  #           bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
  #         }
  #       }
  #       ###message(i," bParNeighbors: ",bParentsOrNeighbors)
  #       ####Then check A and C are not adjacent and B\in S_AC
  #       if(length(bAdjs)>0 && length(bParentsOrNeighbors)>0){
  #         for (k in 1:length(bParentsOrNeighbors)) {
  #           for(l in 1:length(bAdjs)){
  #             aIndex<-bParentsOrNeighbors[k]
  #             cIndex<-bAdjs[l]
  #             ##if A,C are not adjacents and B\in S_AC
  #             if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && is.element(vset[i],sepset(aIndex,cIndex))){
  #               # ##message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
  #               # ##print(skel@sepset[[aIndex]][[cIndex]])
  #               # ##print(skel@sepset[[cIndex]][[aIndex]])
  #               ###then block from C to B
  #               wmat[cIndex,i]<-5
  #               ###message("b: ",vset[i]," ","neighbor: ",vset[nbIndex])
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 2
  #     ########################################################################
  #     ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
  #     ########################################################################
  #     if(lcc>0){
  #       for (i in 1:lcc) {
  #         cc<-ccycles[[i]]
  #         ###print(cc)
  #         #cc[j] plays the role of A
  #         for (j in 1:length(cc)) {
  #           aAdjs<-c()
  #           for (k in 1:length(cc)) {
  #             #aAdjs<-c()
  #             if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
  #               aAdjs<-c(aAdjs,cc[k])
  #             }
  #           }
  #           ###message(cc[j]," as A: ",aAdjs)
  #           if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
  #             #aAdjs[1] is a candidate for being B
  #             bPN<-findPNs(aAdjs[1],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[1],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #           if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
  #             #aAdjs[2] is a candidate for being B
  #             bPN<-findPNs(aAdjs[2],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[2],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 3
  #     
  #     ########################################################################
  #     ###### Rule 4 ##########################################################
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       aAdjs<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           aAdjs<-c(aAdjs,j)
  #         }
  #       }
  #       if(length(aAdjs)>=3){
  #         for (k in 1:length(aAdjs)) {
  #           for (l in 1:length(aAdjs)) {
  #             cIndex<-aAdjs[k]
  #             dIndex<-aAdjs[l]
  #             ##if C,D are not adjacents and A\in S_CD
  #             if((cIndex!=dIndex) && is.element(vset[i],sepset(dIndex,cIndex))){
  #               for (r in 1:length(aAdjs)) {
  #                 bIndex<-aAdjs[r]
  #                 #and if B to C,D is blocked and B to A is not blocked
  #                 if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
  #                   #then block B to A i.e., A +- B
  #                   wmat[bIndex,i]<-5
  #                   change<-TRUE
  #                   break
  #                 }
  #               }
  #             }
  #           }
  #         }
  #       }
  #     }##end of rule 4
  #   } 
  ##################################################################################################
  ################## Replace every edge A +- B with A -> B #########################################
  ################## Replace every edge A +-+ B with A - B #########################################
  ##################################################################################################
  ###print(wmat)
  #wmat<-t(wmat)
  for (i in 1:(length(V)-1)) {
    for (j in (i+1):length(V)) {
      if(wmat[i,j]==5 && wmat[j,i]==5){
        wmat[i,j]<-wmat[j,i]<-1
      }
      if(wmat[i,j]==5 && wmat[j,i]==1){
        wmat[i,j]<-0
        #wmat[j,i]<-0
      }
      if(wmat[j,i]==5 && wmat[i,j]==1){
        wmat[j,i]<-0
        #wmat[i,j]<-0
      }
    }
  }
  return(wmat)
}
####################################################################################################################
`learn.amp.normPC`<-function(data,p.value,method ="stable",...){
  #### First, load the package pcalg and the data set. ####
  V <- colnames(data)
  #covariance<-cov(data)
  #print(covariance)
  n<-nrow(data)
  suffStat<-list(C=cor(data),n=n)
  mm <- switch(method,
               stable = 1,
               original = 2,
               stable.fast = 3,
               0)
  if (mm == 0) stop("Invalid method!")
  if(mm==1){
    skel<-Gskeleton(suffStat,
                    indepTest = gaussCItest, ## (partial correlations)
                    labels =V,alpha=p.value, verbose = FALSE,...)
  }
  if(mm==2){
    skel<-Gskeleton(suffStat,labels =V,method="original",
                    indepTest = gaussCItest, ## (partial correlations)
                    alpha=p.value, verbose = FALSE,...)
  }
  if(mm==3){
    skel<-Gskeleton(suffStat,labels =V,method="stable.fast",
                    indepTest = gaussCItest, ## (partial correlations)
                    alpha=p.value, verbose = FALSE,...)
  }
  wmat<-as(skel@graph,"matrix")####gives adjacency matrix of pcalg, which is corresponding to the learned skeleton
  #print(skel@sepset)
  
  #################################################################
  #################################################################
  #finding chordless cycles 
  ccycles<-ChordlessCycles(wmat)
  lcc<-length(ccycles)
  ########################################################################
  ########## Finding possible parent-neighbors nodes for a given node in 
  ########## a chordless cycle ###########################################
  ########################################################################
  `findPNs`<-function(v,ccycle){
    vParentsNeighbors<-c()
    for (i in 1:length(ccycle)) {
      if(wmat[v,ccycle[i]]==5){
        vParentsNeighbors<-c(vParentsNeighbors,ccycle[i])
      }
    }
    return(vParentsNeighbors)
  }
  ##########################################################################################################
  ##########################################################################################################
  ###A function that determines the type of chordless cycles: with or without blocks########################
  ### if the cycle is chordless and without blocks, the function returns 1 o.w., returns 0 #################
  ##########################################################################################################
  ##########################################################################################################
  ccycles_type<-function(vect,amat){
    cu<-cd<-co<-0
    l<-length(vect)
    for (i in 1:(l-1)) {
      if(amat[vect[i],vect[i+1]]==1 &&  amat[vect[i+1],vect[i]]==1){
        cu<-cu+1
      }
    }
    if(amat[vect[l],vect[1]]==1 &&  amat[vect[1],vect[l]]==1){
      cu<-cu+1
    }
    if(cu==l){
      #undirected cycle
      return(1)
    }else{
      return(0)
    }
  }
  ########################################################################################
  ########################################################################################
  ########## Then apply the rules: 1-4 ###################################################
  ########################################################################################
  change<-TRUE
  while(change){
    change<-FALSE
    ########################################################################
    ###### Rule 1 ###### A o--o B o--o C & B\not\in S_AC -> A +-o B o-+ C
    ########################################################################
    for (i in 1:length(V)) {
      #message(V[i],": *")
      bAdjs<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          bAdjs<-c(bAdjs,j)
        }
      }
      ####Then check A and C are not adjacent and B\not\in S_AC
      if(length(bAdjs)>1){
        for (k in 1:(length(bAdjs)-1)) {
          for(l in (k+1):length(bAdjs)){
            aIndex<-bAdjs[k]
            cIndex<-bAdjs[l]
            # message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
            # print(skel@sepset[[aIndex]][[cIndex]])
            # print(skel@sepset[[cIndex]][[aIndex]])
            ##if they are not adjacents and B\not\in S_AC
            if((wmat[i,aIndex]!=5 || wmat[i,cIndex]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && !((i %in% skel@sepset[[aIndex]][[cIndex]]) ||
                                                                                                                 (i %in% skel@sepset[[cIndex]][[aIndex]]))){
              ###then block from B to A and C
              wmat[i,aIndex]<-5
              wmat[i,cIndex]<-5
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 1
    ########################################################################
    ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
    ########################################################################
    for (i in 1:length(V)) {
      #message("change= ",change)
      bAdjs<-c()#parent, neighbor, or children of B
      bParentsOrNeighbors<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          bAdjs<-c(bAdjs,j)
        }
        if(wmat[i,j]==5){
          bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
        }
      }
      #message(i," bParNeighbors: ",bParentsOrNeighbors)
      ####Then check A and C are not adjacent and B\in S_AC
      if(length(bAdjs)>0 && length(bParentsOrNeighbors)>0){
        for (k in 1:length(bParentsOrNeighbors)) {
          for(l in 1:length(bAdjs)){
            aIndex<-bParentsOrNeighbors[k]
            cIndex<-bAdjs[l]
            ##if A,C are not adjacents and B\in S_AC
            if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && ((i %in% skel@sepset[[aIndex]][[cIndex]]) ||
                                                                                                               (i %in% skel@sepset[[cIndex]][[aIndex]]))){
              # message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
              # print(skel@sepset[[aIndex]][[cIndex]])
              # print(skel@sepset[[cIndex]][[aIndex]])
              ###then block from C to B
              wmat[cIndex,i]<-5
              #message("b: ",V[i]," ","neighbor: ",V[nbIndex])
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 2
    ########################################################################
    ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
    ########################################################################
    if(lcc>0){
      for (i in 1:lcc) {
        cc<-ccycles[[i]]
        #print(cc)
        #cc[j] plays the role of A
        for (j in 1:length(cc)) {
          aAdjs<-c()
          for (k in 1:length(cc)) {
            #aAdjs<-c()
            if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
              aAdjs<-c(aAdjs,cc[k])
            }
          }
          #message(cc[j]," as A: ",aAdjs)
          if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
            #aAdjs[1] is a candidate for being B
            bPN<-findPNs(aAdjs[1],cc)
            if(length(bPN)==1){
              wmat[aAdjs[1],cc[j]]<-5
              change<-TRUE
              break
            }
          }
          if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
            #aAdjs[2] is a candidate for being B
            bPN<-findPNs(aAdjs[2],cc)
            if(length(bPN)==1){
              wmat[aAdjs[2],cc[j]]<-5
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 3
    
    ########################################################################
    ###### Rule 4 ##########################################################
    ########################################################################
    for (i in 1:length(V)) {
      aAdjs<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          aAdjs<-c(aAdjs,j)
        }
      }
      if(length(aAdjs)>=3){
        for (k in 1:length(aAdjs)) {
          for (l in 1:length(aAdjs)) {
            cIndex<-aAdjs[k]
            dIndex<-aAdjs[l]
            ##if C,D are not adjacents and A\in S_CD
            if((cIndex!=dIndex) && ((i %in% skel@sepset[[cIndex]][[dIndex]]) ||
                                    (i %in% skel@sepset[[dIndex]][[cIndex]]))){
              for (r in 1:length(aAdjs)) {
                bIndex<-aAdjs[r]
                #and if B to C,D is blocked and B to A is not blocked
                if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
                  #then block B to A i.e., A +- B
                  wmat[bIndex,i]<-5
                  change<-TRUE
                  break
                }
              }
            }
          }
        }
      }
    }##end of rule 4
  }
  # ##################################################################################################
  #################################################################################################### 
  # ################### This part of code can be used to obtain the Essential graph ##################
  # ################### of the Markov equivalence class of G and not just a CG in class ##############
  # ##### Learning marginal AMP chain graphs under faithfulness revisited (Pena & Gomez-Olmedo, 2016)#
  # ##################################################################################################
  # #  Replace every edge A - B in every cycle that is of length greater than three, chordless,####### 
  # #  and without blocks with A +-+ B ###############################################################
  # ##################################################################################################
  #   if(lcc>0){
  #     for (i in 1:lcc) {
  #       c<-ccycles[i]
  #       l<-length(c)
  #       if((length(c)>3) && (ccycles_type(c,wmat)==1)){
  #         for (j in 1:(l-1)) {
  #           wmat[c[j],c[j+1]]<-wmat[c[j+1],c[j]]<-5
  #         }
  #         wmat[c[1],c[l]]<-wmat[c[l],c[1]]<-5
  #       }
  #     }
  #   }
  # ##################################################################################################
  # ##################################################################################################
  # ##  Apply the rules R2-R4 while possible #########################################################
  # ##################################################################################################
  # ##################################################################################################
  #   change<-TRUE
  #   while(change){
  #     change<-FALSE
  #     ########################################################################
  #     ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       #message("change= ",change)
  #       bAdjs<-c()#parent, neighbor, or children of B
  #       bParentsOrNeighbors<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           bAdjs<-c(bAdjs,j)
  #         }
  #         if(wmat[i,j]==5){
  #           bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
  #         }
  #       }
  #       #message(i," bParNeighbors: ",bParentsOrNeighbors)
  #       ####Then check A and C are not adjacent and B\in S_AC
  #       if(length(bAdjs)>0 && length(bParentsOrNeighbors)>0){
  #         for (k in 1:length(bParentsOrNeighbors)) {
  #           for(l in 1:length(bAdjs)){
  #             aIndex<-bParentsOrNeighbors[k]
  #             cIndex<-bAdjs[l]
  #             ##if A,C are not adjacents and B\in S_AC
  #             if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && ((i %in% skel@sepset[[aIndex]][[cIndex]]) ||
  #                                                                                                                (i %in% skel@sepset[[cIndex]][[aIndex]]))){
  #               # message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
  #               # print(skel@sepset[[aIndex]][[cIndex]])
  #               # print(skel@sepset[[cIndex]][[aIndex]])
  #               ###then block from C to B
  #               wmat[cIndex,i]<-5
  #               #message("b: ",V[i]," ","neighbor: ",V[nbIndex])
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 2
  #     ########################################################################
  #     ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
  #     ########################################################################
  #     if(lcc>0){
  #       for (i in 1:lcc) {
  #         cc<-ccycles[[i]]
  #         #print(cc)
  #         #cc[j] plays the role of A
  #         for (j in 1:length(cc)) {
  #           aAdjs<-c()
  #           for (k in 1:length(cc)) {
  #             #aAdjs<-c()
  #             if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
  #               aAdjs<-c(aAdjs,cc[k])
  #             }
  #           }
  #           #message(cc[j]," as A: ",aAdjs)
  #           if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
  #             #aAdjs[1] is a candidate for being B
  #             bPN<-findPNs(aAdjs[1],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[1],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #           if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
  #             #aAdjs[2] is a candidate for being B
  #             bPN<-findPNs(aAdjs[2],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[2],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 3
  #     
  #     ########################################################################
  #     ###### Rule 4 ##########################################################
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       aAdjs<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           aAdjs<-c(aAdjs,j)
  #         }
  #       }
  #       if(length(aAdjs)>=3){
  #         for (k in 1:length(aAdjs)) {
  #           for (l in 1:length(aAdjs)) {
  #             cIndex<-aAdjs[k]
  #             dIndex<-aAdjs[l]
  #             ##if C,D are not adjacents and A\in S_CD
  #             if((cIndex!=dIndex) && ((i %in% skel@sepset[[cIndex]][[dIndex]]) ||
  #                                     (i %in% skel@sepset[[dIndex]][[cIndex]]))){
  #               for (r in 1:length(aAdjs)) {
  #                 bIndex<-aAdjs[r]
  #                 #and if B to C,D is blocked and B to A is not blocked
  #                 if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
  #                   #then block B to A i.e., A +- B
  #                   wmat[bIndex,i]<-5
  #                   change<-TRUE
  #                   break
  #                 }
  #               }
  #             }
  #           }
  #         }
  #       }
  #     }##end of rule 4
  #   } 
  ##################################################################################################
  ################## Replace every edge A +- B with A -> B #########################################
  ################## Replace every edge A +-+ B with A - B #########################################
  ##################################################################################################
  for (i in 1:(length(V)-1)) {
    for (j in (i+1):length(V)) {
      if(wmat[i,j]==5 && wmat[j,i]==5){
        wmat[i,j]<-wmat[j,i]<-1
      }
      if(wmat[i,j]==5 && wmat[j,i]==1){
        wmat[i,j]<-0
      }
      if(wmat[j,i]==5 && wmat[i,j]==1){
        wmat[j,i]<-0
      }
    }
  }
  return(wmat)
}
#################################################################################################################
## adapt code from the "pcAlgo" from "pcalg" package in R
`.get.exed.cand1` <- function(tree, amat){
  vset <- rownames(amat)
  n.clique <- length(tree@cliques)
  cand.pairs <- c()
  if(n.clique == 1) return(cand.pairs)
  for (i in 1:(n.clique - 1)) {
    sepset <- tree@separators[[i]]$separator
    idx <- match(sepset, vset)
    G <- amat[idx, idx]
    ind <- which(G == 1, arr.ind = TRUE)
    if (any(ind)) {
      ind[,1] <- idx[ind[,1]]
      ind[,2] <- idx[ind[,2]]
      cand.pairs <- rbind(cand.pairs, ind)
    }
  }
  unique(cand.pairs)
}
#######################################################################
`.getNextSet` <- function (n, k, set) 
{
  chInd <- k - (zeros <- sum((seq(n - k + 1, n) - set) == 0))
  wasLast <- (chInd == 0)
  if (!wasLast) {
    set[chInd] <- set[chInd] + 1
    if (chInd < k) 
      set[(chInd + 1):k] <- seq(set[chInd] + 1, set[chInd] + 
                                  zeros)
  }
  list(nextSet = set, wasLast = wasLast)
}
######################################################################

`.get.localug.opcAMPCGs` <- function(cov, n, p.value){
  p <- nrow(as.matrix(cov))
  vset <- rownames(cov)
  sep.pairs <- c()
  G <- matrix(rep(TRUE, p*p), p, p)
  diag(G) <- FALSE
  seq_p <- 1:p
  done <- FALSE
  ord <- 0
  while (!done && any(G)) {
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ind <- ind[order(ind[ ,1]), ]
    remainingEdgeTests <- nrow(ind)
    #G.l <- split(G, gl(p,p))
    for(i in 1:remainingEdgeTests) {
      x <- ind[i,1]
      y <- ind[i,2]
      if (G[y, x]) {
        nbrsBool <- G[,x]
        nbrsBool[y] <- FALSE
        ###find adjacent variables
        nbrs <- seq_p[nbrsBool]
        copyOfnbrs<-nbrs
        ###message("nbrs= ")
        ###print(nbrs)
        LNbs<-length(nbrs)
        ###also find adjacents of adjacents
        if(LNbs>0){
          copyOfnbrs<-nbrs
          for (k in 1:LNbs) {
            tempNbsBool<-G[,nbrs[k]]
            tempNbsBool[x]<-tempNbsBool[y]<-FALSE
            copyOfnbrs<-union(copyOfnbrs,seq_p[tempNbsBool])
          }
          #nbrs<-sort(copyOfnbrs)
          nbrs<-copyOfnbrs
        }
        ###message("nbrs= ")
        ###print(nbrs)
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord)
            done <- FALSE
          S <- seq(length = ord)
          repeat {
            p.val <- norm.ci.test(cov, n, vset[x], vset[y],
                                  vset[nbrs[S]])$p.value
            # nindtest<-nindtest+1
            # message("nindtest= ",nindtest)
            if (p.val > p.value) {
              G[x, y] <- G[y, x] <- FALSE
              pair <- new("sep.pair", u = vset[x],
                          v = vset[y], s = vset[nbrs[S]])
              sep.pairs <- append(sep.pairs, pair)
              break
            }
            else {
              nextSet <- .getNextSet(length_nbrs, ord, S)
              if (nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1
  }
  sep.pairs
}
#############################################################################################
`learn.original.skeleton.normAMPCGs` <- function(tree, cov, n, p.value, drop = TRUE)
{
  #nindtest<-0
  validObject(tree)
  local.ug <- c()
  vset <- rownames(cov)
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    #        if (length(idx) >= 10)
    new.ug <- .get.localug.opcAMPCGs(cov[idx, idx], n, p.value)
    #        else
    #            new.ug <- .get.localug.ic(cov[idx, idx], n, p.value)
    local.ug <- append(local.ug, new.ug)
  }
  p <- length(vset)
  amat <- matrix(0, p, p)
  rownames(amat) <- colnames(amat) <- vset
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    amat[idx, idx] <- 1
  }
  diag(amat) <- 0
  sep.pairs <- c()
  n.loc.sep <- length(local.ug)
  if(n.loc.sep>0)
    for(i in 1:n.loc.sep){
      u <- local.ug[[i]]@u
      v <- local.ug[[i]]@v
      if(amat[u,v] == 1){
        amat[u,v] <- amat[v,u] <- 0
        sep.pairs <- append(sep.pairs, local.ug[i])
      }
    }
  
  ## the following code is partially adapted from the "pcAlgo" function
  ## from "pcalg" package in R
  
  if (drop) {
    ind <- .get.exed.cand1(tree, amat)
    if (any(ind)) {
      ind <- ind[order(ind[,1]),]
      ord <- 0
      seq_p <- 1:p
      done <- FALSE
      remainingEdgeTests <- nrow(ind)
      while (!done && any(as.logical(amat))) {
        done <- TRUE
        #amat.l <- split(amat, gl(p,p))
        for (i in 1:remainingEdgeTests) {
          x <- ind[i, 1]
          y <- ind[i, 2]
          if (amat[y, x]) {
            nbrsBool <- amat[,x] == 1
            nbrsBool[y] <- FALSE
            ###find adjacent variables
            nbrs <- seq_p[nbrsBool]
            copyOfnbrs<-nbrs
            ###message("nbrs= ")
            ###print(nbrs)
            LNbs<-length(nbrs)
            ###also find adjacents of adjacents
            if(LNbs>0){
              copyOfnbrs<-nbrs
              for (k in 1:LNbs) {
                tempNbsBool<-amat[,nbrs[k]]==1
                tempNbsBool[x]<-tempNbsBool[y]<-FALSE
                copyOfnbrs<-union(copyOfnbrs,seq_p[tempNbsBool])
              }
              #nbrs<-sort(copyOfnbrs)
              nbrs<-copyOfnbrs
            }
            ###message("nbrs= ")
            ###print(nbrs)
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord)
                done <- FALSE
              S <- seq(length = ord)
              repeat {
                p.val <- norm.ci.test(cov, n, vset[x], vset[y],
                                      vset[nbrs[S]])$p.value
                #nindtest<-nindtest+1
                if (p.val > p.value) {
                  amat[x, y] <- amat[y, x] <- 0
                  pair <- new("sep.pair", u = vset[x],
                              v = vset[y], s = vset[nbrs[S]])
                  sep.pairs <- append(sep.pairs, pair)
                  break
                }
                else {
                  nextSet <- .getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          }
        }
        ord <- ord + 1
      }
    }
    ##         } else {
    ##             if (any(ind)) {
    ##                 for(i in 1:nrow(ind)){
    ##                     pair <- new("sep.pair", u = vset[ind[i,1]],
    ##                                 v = vset[ind[i,2]], s = character(0))
    ##                     cand <- setdiff(vset[amat[pair@u,]==1], pair@v)
    ##                     idx <- c(pair@u, pair@v, cand)
    ##                     res <- .get.sep(cov[idx, idx], n, p.value, pair@u, pair@v, cand)
    ##                     if(res$seped){
    ##                         amat[pair@u, pair@v] <- amat[pair@v, pair@u] <- 0
    ##                         sep.pairs <- append(sep.pairs, res$sep)
    ##                     } 
    ##                 }
    ##             }
    ##         }
  }
  #
  #return(nindtest)
  return(list(amat=amat, sep.pairs=sep.pairs))    
}

###########################################################################################################3
######################################################################

`.get.localug.pcAMPCGs` <- function(cov, n, p.value){
  p <- nrow(as.matrix(cov))
  vset <- rownames(cov)
  sep.pairs <- c()
  G <- matrix(rep(TRUE, p*p), p, p)
  diag(G) <- FALSE
  seq_p <- 1:p
  done <- FALSE
  ord <- 0
  while (!done && any(G)) {
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ind <- ind[order(ind[ ,1]), ]
    remainingEdgeTests <- nrow(ind)
    G.l <- split(G, gl(p,p))
    for(i in 1:remainingEdgeTests) {
      x <- ind[i,1]
      y <- ind[i,2]
      if (G[y, x]) {
        nbrsBool <- G.l[[x]]
        nbrsBool[y] <- FALSE
        ###find adjacent variables
        nbrs <- seq_p[nbrsBool]
        copyOfnbrs<-nbrs
        ###message("nbrs= ")
        ###print(nbrs)
        LNbs<-length(nbrs)
        ###also find adjacents of adjacents
        if(LNbs>0){
          copyOfnbrs<-nbrs
          for (k in 1:LNbs) {
            tempNbsBool<-G.l[[nbrs[k]]]
            tempNbsBool[x]<-tempNbsBool[y]<-FALSE
            copyOfnbrs<-union(copyOfnbrs,seq_p[tempNbsBool])
          }
          #nbrs<-sort(copyOfnbrs)
          nbrs<-copyOfnbrs
        }
        ###message("nbrs= ")
        ###print(nbrs)
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord)
            done <- FALSE
          S <- seq(length = ord)
          repeat {
            p.val <- norm.ci.test(cov, n, vset[x], vset[y],
                                  vset[nbrs[S]])$p.value
            #nindtest<-nindtest+1
            if (p.val > p.value) {
              G[x, y] <- G[y, x] <- FALSE
              pair <- new("sep.pair", u = vset[x],
                          v = vset[y], s = vset[nbrs[S]])
              sep.pairs <- append(sep.pairs, pair)
              break
            }
            else {
              nextSet <- .getNextSet(length_nbrs, ord, S)
              if (nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1
  }
  sep.pairs
}
#############################################################################################
`learn.stable.skeleton.normAMPCGs` <- function(tree, cov, n, p.value, drop = TRUE)
{
  validObject(tree)
  local.ug <- c()
  vset <- rownames(cov)
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    #        if (length(idx) >= 10)
    new.ug <- .get.localug.pcAMPCGs(cov[idx, idx], n, p.value)
    #        else
    #            new.ug <- .get.localug.ic(cov[idx, idx], n, p.value)
    local.ug <- append(local.ug, new.ug)
  }
  p <- length(vset)
  amat <- matrix(0, p, p)
  rownames(amat) <- colnames(amat) <- vset
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    amat[idx, idx] <- 1
  }
  diag(amat) <- 0
  sep.pairs <- c()
  n.loc.sep <- length(local.ug)
  if(n.loc.sep>0)
    for(i in 1:n.loc.sep){
      u <- local.ug[[i]]@u
      v <- local.ug[[i]]@v
      if(amat[u,v] == 1){
        amat[u,v] <- amat[v,u] <- 0
        sep.pairs <- append(sep.pairs, local.ug[i])
      }
    }
  
  ## the following code is partially adapted from the "pcAlgo" function
  ## from "pcalg" package in R
  
  if (drop) {
    ind <- .get.exed.cand1(tree, amat)
    if (any(ind)) {
      ind <- ind[order(ind[,1]),]
      ord <- 0
      seq_p <- 1:p
      done <- FALSE
      remainingEdgeTests <- nrow(ind)
      while (!done && any(as.logical(amat))) {
        done <- TRUE
        amat.l <- split(amat, gl(p,p))
        for (i in 1:remainingEdgeTests) {
          x <- ind[i, 1]
          y <- ind[i, 2]
          if (amat[y, x]) {
            nbrsBool <- amat.l[[x]] == 1
            nbrsBool[y] <- FALSE
            ###find adjacent variables
            nbrs <- seq_p[nbrsBool]
            copyOfnbrs<-nbrs
            ###message("nbrs= ")
            ###print(nbrs)
            LNbs<-length(nbrs)
            ###also find adjacents of adjacents
            if(LNbs>0){
              copyOfnbrs<-nbrs
              for (k in 1:LNbs) {
                tempNbsBool<-amat.l[[nbrs[k]]]==1
                tempNbsBool[x]<-tempNbsBool[y]<-FALSE
                copyOfnbrs<-union(copyOfnbrs,seq_p[tempNbsBool])
              }
              #nbrs<-sort(copyOfnbrs)
              nbrs<-copyOfnbrs
            }
            ###message("nbrs= ")
            ###print(nbrs)
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord)
                done <- FALSE
              S <- seq(length = ord)
              repeat {
                p.val <- norm.ci.test(cov, n, vset[x], vset[y],
                                      vset[nbrs[S]])$p.value
                #nindtest<-nindtest+1
                if (p.val > p.value) {
                  amat[x, y] <- amat[y, x] <- 0
                  pair <- new("sep.pair", u = vset[x],
                              v = vset[y], s = vset[nbrs[S]])
                  sep.pairs <- append(sep.pairs, pair)
                  break
                }
                else {
                  nextSet <- .getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          }
        }
        ord <- ord + 1
      }
    }
    ##         } else {
    ##             if (any(ind)) {
    ##                 for(i in 1:nrow(ind)){
    ##                     pair <- new("sep.pair", u = vset[ind[i,1]],
    ##                                 v = vset[ind[i,2]], s = character(0))
    ##                     cand <- setdiff(vset[amat[pair@u,]==1], pair@v)
    ##                     idx <- c(pair@u, pair@v, cand)
    ##                     res <- .get.sep(cov[idx, idx], n, p.value, pair@u, pair@v, cand)
    ##                     if(res$seped){
    ##                         amat[pair@u, pair@v] <- amat[pair@v, pair@u] <- 0
    ##                         sep.pairs <- append(sep.pairs, res$sep)
    ##                     } 
    ##                 }
    ##             }
    ##         }
  }
  #
  #return(nindtest)
  return(list(amat=amat, sep.pairs=sep.pairs))    
}

#############################################################################################################
######################################################################
.get.localug.m.ic <- function(freq.tb, p.value)
{
  p <- length(freq.tb@levels)
  vset <- colnames(freq.tb@table)
  vset <- vset[-length(vset)]
  sep.pairs <- c()
  if (p > 1)
    for(i in 1:(p-1))
      for(j in (i+1):p){
        cand <- vset[-c(i,j)]
        res <- .get.sep.m(freq.tb, p.value, vset[i], vset[j], cand)
        if(res$seped)
          sep.pairs <- append(sep.pairs, res$sep)
      }
  sep.pairs
}

################################################################################
.get.localug.m.pcAMPCGs <- function(freq.tb, p.value)
{
  p <- length(freq.tb@levels)
  vset <- colnames(freq.tb@table)
  vset <- vset[-length(vset)]
  sep.pairs <- c()
  G <- matrix(rep(TRUE, p*p), p, p)
  diag(G) <- FALSE
  seq_p <- 1:p
  done <- FALSE
  ord <- 0
  n.edgetests <- numeric(1)
  while (!done && any(G)) {
    n.edgetests[ord + 1] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ind <- ind[order(ind[ ,1]), ]
    remainingEdgeTests <- nrow(ind)
    G.l <- split(G, gl(p,p))
    for(i in 1:remainingEdgeTests) {
      x <- ind[i,1]
      y <- ind[i,2]
      if (G[y, x]) {
        nbrsBool <- G.l[[x]]
        nbrsBool[y] <- FALSE
        ###find adjacent variables
        nbrs <- seq_p[nbrsBool]
        copyOfnbrs<-nbrs
        ###message("nbrs= ")
        ###print(nbrs)
        LNbs<-length(nbrs)
        ###also find adjacents of adjacents
        if(LNbs>0){
          copyOfnbrs<-nbrs
          for (k in 1:LNbs) {
            tempNbsBool<-G.l[[nbrs[k]]]
            tempNbsBool[x]<-tempNbsBool[y]<-FALSE
            copyOfnbrs<-union(copyOfnbrs,seq_p[tempNbsBool])
          }
          #nbrs<-sort(copyOfnbrs)
          nbrs<-copyOfnbrs
        }
        ###message("nbrs= ")
        ###print(nbrs)
        length_nbrs <- length(nbrs)
        if (length_nbrs >= ord) {
          if (length_nbrs > ord)
            done <- FALSE
          S <- seq(length = ord)
          repeat {
            n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
            p.val <- multinom.ci.test(freq.tb, vset[x], vset[y],
                                      vset[nbrs[S]])$p.value
            if (p.val > p.value) {
              G[x, y] <- G[y, x] <- FALSE
              pair <- new("sep.pair", u = vset[x],
                          v = vset[y], s = vset[nbrs[S]])
              sep.pairs <- append(sep.pairs, pair)
              break
            }
            else {
              nextSet <- .getNextSet(length_nbrs, ord, S)
              if (nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          }
        }
      }
    }
    ord <- ord + 1
  }
  sep.pairs
}

##############################################################################################

.get.sep.m <- function(freq.tb, p.value, u, v, cand)
{
  seped <- FALSE
  p <- length(cand)
  if (p < 8) {
    mat <- as.matrix(t(sapply(1:(2^p),
                              function(x) .binary(x-1, max(1,p))$dicotomy)))
    if(nrow(mat) < ncol(mat)) mat <- t(mat)
    p.val <- apply(mat, 1, function(x) multinom.ci.test(freq.tb, u, v, cand[x])$p.value)
    p.val.max <- max(p.val)
    idx <- which(p.val == p.val.max)
    sep <- new("sep.pair",
               u=u,
               v=v,
               s=character(0))
    if(p.val.max >= p.value){
      sep@s <- cand[mat[idx[1],]]
      seped <- TRUE
    }
    return(list(seped=seped, sep=sep))
  } else {
    .get.sep.m.stepwise(freq.tb, p.value, u, v, cand)
  }
}


.get.sep.m.stepwise <- function(freq.tb, p.value, u, v, cand){
  .get.sep.m.step(freq.tb, p.value, u, v, cand, c())
}

.get.sep.m.step <- function(freq.tb, p.value, u, v, current, rest)
{
  modified <- FALSE
  sep <- new("sep.pair",
             u = u,
             v = v,
             s = character(0))
  pp <- multinom.ci.test(freq.tb, u, v, current)$p.value
  if(pp >= p.value) {
    sep@s = current
    return(list(seped = TRUE, sep = sep))
  }
  if(length(current) == 0)
    return(list(seped = FALSE, sep = sep))
  n.curr <- length(current)
  n.rest <- length(rest)
  mat <- matrix(rep(current, n.curr), n.curr, byrow = TRUE)
  diag(mat) <- NA
  pval <- apply(mat, 1, function(x) multinom.ci.test(freq.tb, u, v, x[!is.na(x)])$p.value)
  todel <- which(pval == max(pval))
  if(length(todel) > 1)
    todel <- todel[sample(length(todel),1)]
  if(pval[todel] >= pp){
    todel <- current[todel]
    current <- setdiff(current, todel)
    rest <- union(rest, todel)
    pp <- max(pval)
    modified <- TRUE
  }
  if(pp >= p.value){
    sep@s = current
    return(list(seped = TRUE, sep = sep))
  }
  if(length(rest) == 0)
    return(list(seped = FALSE, sep = sep))
  if(modified){
    n.curr <- length(current)
    n.rest <- length(rest)
  }
  mat <- matrix(NA, n.rest, n.curr + 1)
  if(n.curr != 0)
    mat[,1:n.curr] <- matrix(rep(current, n.rest), n.rest, byrow = TRUE)
  mat[,n.curr + 1] <- rest
  pval <- apply(mat, 1, function(x) multinom.ci.test(freq.tb, u, v, x)$p.value)
  toadd <- which(pval == max(pval))
  if(length(toadd) > 1)
    toadd <- toadd[sample(length(toadd),1)]
  if(pval[toadd] > pp){
    toadd <- rest[toadd]
    current <- union(current, toadd)
    rest <- setdiff(rest, toadd)
    modified <- TRUE
  }
  if(modified)
    return(.get.sep.m.step(freq.tb, p.value, u, v, current, rest))
  else 
    return(list(seped = FALSE, sep = sep))
}
########################################################################################

`learn.skeleton.multinomAMPCGs` <- function(tree, freq.tb, p.value, drop = TRUE)
{
  validObject(tree)
  validObject(freq.tb)
  local.ug <- c()
  vset <- colnames(freq.tb@table)
  vset <- vset[-length(vset)]
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    #        if (length(idx) >= 10)
    new.ug <- .get.localug.m.pcAMPCGs(compress.freq.tb(freq.tb, idx), p.value)
    #        else
    #            new.ug <- .get.localug.m.ic(compress.freq.tb(freq.tb, idx), p.value)
    local.ug <- append(local.ug, new.ug)
  }
  p <- length(vset)
  amat <- matrix(0, p, p)
  rownames(amat) <- colnames(amat) <- vset
  n.clique <- length(tree@cliques)
  for(i in 1:n.clique){
    idx <- tree@cliques[[i]]$vset
    amat[idx, idx] <- 1
  }
  diag(amat) <- 0
  sep.pairs <- c()
  n.loc.sep <- length(local.ug)
  if(n.loc.sep>0)
    for(i in 1:n.loc.sep){
      u <- local.ug[[i]]@u
      v <- local.ug[[i]]@v
      if(amat[u,v] == 1){
        amat[u,v] <- amat[v,u] <- 0
        sep.pairs <- append(sep.pairs, local.ug[i])
      }
    }
  
  
  ## code partially adapted from the "pcAlgo" function from the
  ## "pcalg" package
  
  if (drop) {
    ind <- .get.exed.cand1(tree, amat)
    #        if (any(ind) && nrow(ind) > 50) {
    if (any(ind)) {
      ind <- ind[order(ind[,1]),]
      ord <- 0
      seq_p <- 1:p
      done <- FALSE
      remainingEdgeTests <- nrow(ind)
      while (!done && any(as.logical(amat))) {
        done <- TRUE
        amat.l <- split(amat, gl(p,p))
        for (i in 1:remainingEdgeTests) {
          x <- ind[i, 1]
          y <- ind[i, 2]
          if (amat[y, x]) {
            nbrsBool <- amat.l[[x]] == 1
            nbrsBool[y] <- FALSE
            ###find adjacent variables
            nbrs <- seq_p[nbrsBool]
            copyOfnbrs<-nbrs
            ###message("nbrs= ")
            ###print(nbrs)
            LNbs<-length(nbrs)
            ###also find adjacents of adjacents
            if(LNbs>0){
              copyOfnbrs<-nbrs
              for (k in 1:LNbs) {
                tempNbsBool<-amat.l[[nbrs[k]]]==1
                tempNbsBool[x]<-tempNbsBool[y]<-FALSE
                copyOfnbrs<-union(copyOfnbrs,seq_p[tempNbsBool])
              }
              #nbrs<-sort(copyOfnbrs)
              nbrs<-copyOfnbrs
            }
            ###message("nbrs= ")
            ###print(nbrs)
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord)
                done <- FALSE
              S <- seq(length = ord)
              repeat {
                p.val <- multinom.ci.test(freq.tb, vset[x], vset[y],
                                          vset[nbrs[S]])$p.value
                if (p.val > p.value) {
                  amat[x, y] <- amat[y, x] <- 0
                  pair <- new("sep.pair", u = vset[x],
                              v = vset[y], s = vset[nbrs[S]])
                  sep.pairs <- append(sep.pairs, pair)
                  break
                }
                else {
                  nextSet <- .getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          }
        }
        ord <- ord + 1
      }
    }
    ##         } else {
    ##             if (any(ind)) {
    ##                 for(i in 1:nrow(ind)){
    ##                     pair <- new("sep.pair", u = vset[ind[i,1]],
    ##                                 v = vset[ind[i,2]], s = character(0))
    ##                     cand <- setdiff(vset[amat[pair@u,]==1], pair@v)
    ##                     idx <- c(pair@u, pair@v, cand)
    ##                     res <- .get.sep.m(compress.freq.tb(freq.tb, idx), p.value, pair@u, pair@v, cand)
    ##                     if(res$seped){
    ##                         amat[pair@u, pair@v] <- amat[pair@v, pair@u] <- 0
    ##                         sep.pairs <- append(sep.pairs, res$sep)
    ##                     } 
    ##                 }
    ##             }
    ##         }
  }
  return(list(amat=amat, sep.pairs=sep.pairs))    
}

###############################################################################################################
`learn.original.amp.normLCD`<-function(tgdata,p.value){
  #### First, load the package pcalg and the data set. ####
  V <- colnames(tgdata)
  tgug <- naive.getug.norm(tgdata, p.value)
  tg.jtree <- ug.to.jtree(tgug)
  skel<-learn.original.skeleton.normAMPCGs(tg.jtree, cov(tgdata), n=nrow(tgdata), p.value, drop = TRUE)
  wmat <- skel$amat#gives adjacency matrix of the learned skeleton via LCD algorithm
  vset <- rownames(wmat)
  sep.pairs <- skel$sep.pairs
  cliques <- tg.jtree@cliques
  n.clique <- length(cliques)
  ###print(wmat)
  ###print(tg.jtree)
  #################################################################
  #################################################################
  #finding chordless cycles 
  ccycles<-ChordlessCycles(wmat)
  lcc<-length(ccycles)
  #finding chordless cycles of length 3
  ccyclesOf3<-list()
  vtemp<-c()
  if(lcc>0){
    for (i in 1:lcc){
      if(length(ccycles[[i]])==3){
        vtemp<-c(vtemp,i)
      }
    }
  }
  if(length(vtemp)>0){
    ccyclesOf3<-ccycles[vtemp]
  }
  lccOf3<-length(ccyclesOf3)
  ########################################################################
  ########## Finding possible parent-neighbors nodes for a given node in 
  ########## a chordless cycle ###########################################
  ########################################################################
  `findPNs`<-function(v,ccycle){
    vParentsNeighbors<-c()
    for (i in 1:length(ccycle)) {
      if(wmat[v,ccycle[i]]==5){
        vParentsNeighbors<-c(vParentsNeighbors,ccycle[i])
      }
    }
    return(vParentsNeighbors)
  }
  # ##########################################################################################################
  # ##########################################################################################################
  # ###A function that determines the type of chordless cycles: with or without blocks########################
  # ### if the cycle is chordless and without blocks, the function returns 1 o.w., returns 0 #################
  # ##########################################################################################################
  # ##########################################################################################################
  # ccycles_type<-function(vect,amat){
  #   cu<-cd<-co<-0
  #   l<-length(vect)
  #   for (i in 1:(l-1)) {
  #     if(amat[vect[i],vect[i+1]]==1 &&  amat[vect[i+1],vect[i]]==1){
  #       cu<-cu+1
  #     }
  #   }
  #   if(amat[vect[l],vect[1]]==1 &&  amat[vect[1],vect[l]]==1){
  #     cu<-cu+1
  #   }
  #   if(cu==l){
  #     #undirected cycle
  #     return(1)
  #   }else{
  #     return(0)
  #   }
  # }
  # ########################################################################################
  # ########## A function that returns the separation set of two given variables ###########
  # ########################################################################################
  `AsepCD`<-function(u,v,w){
    pair <- new("sep.pair", u=vset[u], v=vset[v])
    a <- lapply(sep.pairs, function(x) all.equal(x,pair))
    ###message(a," for ",u," and ", v)
    if(length((which(a==TRUE)))>0){
      sep<-sep.pairs[[which(a==TRUE)[1]]]@s
      if(is.element(vset[w], sep)){
        return(1)
      }else{
        return(0)
      }
    }else{
      return(2)
    }
  }
  
  ########################################################################################
  ########################################################################################
  ########## Then apply the rules: 1-4 ###################################################
  ########################################################################################
  change<-TRUE
  while(change){
    change<-FALSE
    ########################################################################
    ###### Rule 1 ###### A o--o B o--o C & B\not\in S_AC -> A +-o B o-+ C
    ########################################################################
    for(i in 1:n.clique){
      cvset <- cliques[[i]]$vset
      ##print(cvset)
      p <- length(cvset)
      if(p > 2){
        for(j in 1:(p-1)){
          for(k in (j+1):p){
            for(l in 1:p){
              u <- cvset[j]
              v <- cvset[k]
              w <- cvset[l]
              ####Then check A and C are not adjacent and B\not\in S_AC
              if((wmat[u,v]+wmat[v,u] == 0) &&
                 (wmat[w,u] != 5 || wmat[w,v] != 5) &&
                 (wmat[u,w]+wmat[w,u] > 0) && (wmat[w,v]+wmat[v,w] > 0)){
                pair <- new("sep.pair", u=u, v=v)
                a <- lapply(sep.pairs, function(x) all.equal(x,pair))
                ###message(a," for ",u," and ", v)
                sep <- sep.pairs[[which(a==TRUE)[1]]]@s
                ###message(sep," for ",u," and ", v)
                if(!is.element(w, sep)){
                  ###then block from B to A and C
                  wmat[w,u] <-5
                  wmat[w,v] <-5
                  change<-TRUE
                  ##message(w ," is a sudo-collider for ", u ," and ", v," rule 1")
                  ##print(wmat)
                  break
                }
              }#second if
              #count<- count + 1
            }#fourth for loop
          }#third for loop
        }#second for loop
      }#first if
      ###message("counter= ",count)
    }#first for loop
    #end of rule 1
    ########################################################################
    ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
    ########################################################################
    for (i in 1:length(V)) {
      ###message("change= ",change)
      bAdjs<-c()#parent, neighbor, or children of B
      bParentsOrNeighbors<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          bAdjs<-c(bAdjs,j)
        }
        if(wmat[i,j]==5){
          bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
        }
      }
      ###message(i," bParNeighbors: ",bParentsOrNeighbors)
      ####Then check A and C are not adjacent and B\in S_AC
      if(length(bAdjs)>1 && length(bParentsOrNeighbors)>0){
        for (k in 1:length(bParentsOrNeighbors)) {
          for(l in 1:length(bAdjs)){
            aIndex<-bParentsOrNeighbors[k]
            cIndex<-bAdjs[l]
            ##if A,C are not adjacents and B\in S_AC
            if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && (wmat[aIndex,cIndex]+ wmat[cIndex,aIndex]==0) && (AsepCD(aIndex,cIndex,i)==1 ||
                                                                                                            AsepCD(aIndex,cIndex,i)==2)){
              # ##message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
              # ##print(skel@sepset[[aIndex]][[cIndex]])
              # ##print(skel@sepset[[cIndex]][[aIndex]])
              ###then block from C to B
              wmat[cIndex,i]<-5
              ##message("b: ",vset[i]," "," c: ",vset[cIndex]," a: ",vset[aIndex]," rule 2")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 2
    ########################################################################
    ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
    ########################################################################
    if(lcc>0){
      for (i in 1:lcc) {
        cc<-ccycles[[i]]
        ###print(cc)
        #cc[j] plays the role of A
        for (j in 1:length(cc)) {
          aAdjs<-c()
          for (k in 1:length(cc)) {
            #aAdjs<-c()
            if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
              aAdjs<-c(aAdjs,cc[k])
            }
          }
          ###message(cc[j]," as A: ",aAdjs)
          if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
            #aAdjs[1] is a candidate for being B
            bPN<-findPNs(aAdjs[1],cc)
            if(length(bPN)==1){
              wmat[aAdjs[1],cc[j]]<-5
              ##message("aAdjs[1] ",aAdjs[1]," cc[j] ",cc[j], "rule 3a")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
          if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
            #aAdjs[2] is a candidate for being B
            bPN<-findPNs(aAdjs[2],cc)
            if(length(bPN)==1){
              wmat[aAdjs[2],cc[j]]<-5
              ##message("aAdjs[1] ",aAdjs[1]," cc[j] ",cc[j], "rule 3b")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 3
    ###Now, since rule 1, 2, and 3 are faster than rule 4, apply these until no change, then go on to rule 4
    if (change){
      next
    }
    
    ########################################################################
    ###### Rule 4 ##########################################################
    ########################################################################
    
    if(lccOf3>1){
      for (i in 1:((lccOf3)-1)) {
        for (j in (i+1):lccOf3) {
          abIndex<-intersect(ccyclesOf3[[i]],ccyclesOf3[[j]])
          cIndex<-setdiff(ccyclesOf3[[i]],abIndex)
          dIndex<-setdiff(ccyclesOf3[[j]],abIndex)
          if((length(abIndex)==2) && (wmat[cIndex,dIndex]+wmat[dIndex,cIndex]==0)){
            for (turn in 1:2) {
              aIndex <- if(turn == 1) abIndex[1] else abIndex[2]
              bIndex <- if(turn == 1) abIndex[2] else abIndex[1]
              ##if C,D are not adjacents and A\in S_CD and if B to C,D is blocked and B to A is not blocked
              if(wmat[bIndex,aIndex]!=5 && wmat[aIndex,bIndex]!=0 && wmat[bIndex,cIndex]==5 && wmat[bIndex,dIndex]==5 &&
                 (AsepCD(cIndex,dIndex,aIndex)==1 || AsepCD(cIndex,dIndex,aIndex)==2)){
                #then block B to A i.e., A +- B
                wmat[bIndex,aIndex]<-5
                ##message("bIndex ",vset[bIndex]," aIndex ",vset[aIndex]," cIndex ",vset[cIndex]," dIndex ",vset[dIndex])
                ##print(wmat)
                change<-TRUE
                break
              }
            }
          }
        }
      }
    }##end of rule 4
    
    # for (i in 1:length(V)) {
    #   aAdjs<-c()
    #   for (j in 1:length(V)) {
    #     if(wmat[j,i]+wmat[i,j]>0){
    #       aAdjs<-c(aAdjs,j)
    #     }
    #   }
    #   if(length(aAdjs)>=3){
    #     for (k in 1:length(aAdjs)) {
    #       for (l in 1:length(aAdjs)) {
    #         cIndex<-aAdjs[k]
    #         dIndex<-aAdjs[l]
    #         ##if C,D are not adjacents and A\in S_CD
    #         if((cIndex!=dIndex) && !((wmat[i,cIndex]+wmat[cIndex,i]==2 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,dIndex]+wmat[dIndex,i]==2 && wmat[i,cIndex]==5 && wmat[cIndex,i]==1) ||
    #                                 (wmat[i,cIndex]==5 && wmat[cIndex,i]==1 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,cIndex]+wmat[cIndex,i]==10 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,dIndex]+wmat[dIndex,i]==10 && wmat[i,cIndex]==5 && wmat[cIndex,i]==1))){
    #           for (r in 1:length(aAdjs)) {
    #             bIndex<-aAdjs[r]
    #             #and if B to C,D is blocked and B to A is not blocked
    #             if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
    #               #then block B to A i.e., A +- B
    #               wmat[bIndex,i]<-5
    #               ##message("bIndex ",vset[bIndex]," A ",vset[i]," cIndex ",vset[cIndex]," dIndex ",vset[dIndex])
    #               ##print(aAdjs)
    #               ##print(wmat)
    #               change<-TRUE
    #               #break
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }##end of rule 4
  }
  # ##################################################################################################
  #################################################################################################### 
  # ################### This part of code can be used to obtain the Essential graph ##################
  # ################### of the Markov equivalence class of G and not just a CG in class ##############
  # ##### Learning marginal AMP chain graphs under faithfulness revisited (Pena & Gomez-Olmedo, 2016)#
  # ##################################################################################################
  # ##################################################################################################
  # ##################################################################################################
  # #  Replace every edge A - B in every cycle that is of length greater than three, chordless,####### 
  # #  and without blocks with A +-+ B ###############################################################
  # ##################################################################################################
  #   if(lcc>0){
  #     for (i in 1:lcc) {
  #       c<-ccycles[i]
  #       l<-length(c)
  #       if((length(c)>3) && (ccycles_type(c,wmat)==1)){
  #         for (j in 1:(l-1)) {
  #           wmat[c[j],c[j+1]]<-wmat[c[j+1],c[j]]<-5
  #         }
  #         wmat[c[1],c[l]]<-wmat[c[l],c[1]]<-5
  #       }
  #     }
  #   }
  # ##################################################################################################
  # ##################################################################################################
  # ##  Apply the rules R2-R4 while possible #########################################################
  # ##################################################################################################
  # ##################################################################################################
  #   change<-TRUE
  #   while(change){
  #     change<-FALSE
  #     ########################################################################
  #     ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       ###message("change= ",change)
  #       bAdjs<-c()#parent, neighbor, or children of B
  #       bParentsOrNeighbors<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           bAdjs<-c(bAdjs,j)
  #         }
  #         if(wmat[i,j]==5){
  #           bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
  #         }
  #       }
  #       ###message(i," bParNeighbors: ",bParentsOrNeighbors)
  #       ####Then check A and C are not adjacent and B\in S_AC
  #       if(length(bAdjs)>0 && length(bParentsOrNeighbors)>0){
  #         for (k in 1:length(bParentsOrNeighbors)) {
  #           for(l in 1:length(bAdjs)){
  #             aIndex<-bParentsOrNeighbors[k]
  #             cIndex<-bAdjs[l]
  #             ##if A,C are not adjacents and B\in S_AC
  #             if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && is.element(vset[i],sepset(aIndex,cIndex))){
  #               # ##message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
  #               # ##print(skel@sepset[[aIndex]][[cIndex]])
  #               # ##print(skel@sepset[[cIndex]][[aIndex]])
  #               ###then block from C to B
  #               wmat[cIndex,i]<-5
  #               ###message("b: ",vset[i]," ","neighbor: ",vset[nbIndex])
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 2
  #     ########################################################################
  #     ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
  #     ########################################################################
  #     if(lcc>0){
  #       for (i in 1:lcc) {
  #         cc<-ccycles[[i]]
  #         ###print(cc)
  #         #cc[j] plays the role of A
  #         for (j in 1:length(cc)) {
  #           aAdjs<-c()
  #           for (k in 1:length(cc)) {
  #             #aAdjs<-c()
  #             if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
  #               aAdjs<-c(aAdjs,cc[k])
  #             }
  #           }
  #           ###message(cc[j]," as A: ",aAdjs)
  #           if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
  #             #aAdjs[1] is a candidate for being B
  #             bPN<-findPNs(aAdjs[1],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[1],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #           if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
  #             #aAdjs[2] is a candidate for being B
  #             bPN<-findPNs(aAdjs[2],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[2],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 3
  #     
  #     ########################################################################
  #     ###### Rule 4 ##########################################################
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       aAdjs<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           aAdjs<-c(aAdjs,j)
  #         }
  #       }
  #       if(length(aAdjs)>=3){
  #         for (k in 1:length(aAdjs)) {
  #           for (l in 1:length(aAdjs)) {
  #             cIndex<-aAdjs[k]
  #             dIndex<-aAdjs[l]
  #             ##if C,D are not adjacents and A\in S_CD
  #             if((cIndex!=dIndex) && is.element(vset[i],sepset(dIndex,cIndex))){
  #               for (r in 1:length(aAdjs)) {
  #                 bIndex<-aAdjs[r]
  #                 #and if B to C,D is blocked and B to A is not blocked
  #                 if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
  #                   #then block B to A i.e., A +- B
  #                   wmat[bIndex,i]<-5
  #                   change<-TRUE
  #                   break
  #                 }
  #               }
  #             }
  #           }
  #         }
  #       }
  #     }##end of rule 4
  #   } 
  ##################################################################################################
  ################## Replace every edge A +- B with A -> B #########################################
  ################## Replace every edge A +-+ B with A - B #########################################
  ##################################################################################################
  ###print(wmat)
  #wmat<-t(wmat)
  for (i in 1:(length(V)-1)) {
    for (j in (i+1):length(V)) {
      if(wmat[i,j]==5 && wmat[j,i]==5){
        wmat[i,j]<-wmat[j,i]<-1
      }
      if(wmat[i,j]==5 && wmat[j,i]==1){
        wmat[i,j]<-0
        #wmat[j,i]<-0
      }
      if(wmat[j,i]==5 && wmat[i,j]==1){
        wmat[j,i]<-0
        #wmat[i,j]<-0
      }
    }
  }
  return(wmat)
}
############################################################################################################################
`learn.stable.amp.normLCD`<-function(tgdata,p.value){
  #### First, load the package pcalg and the data set. ####
  V <- colnames(tgdata)
  tgug <- naive.getug.norm(tgdata, p.value)
  tg.jtree <- ug.to.jtree(tgug)
  skel<-learn.stable.skeleton.normAMPCGs(tg.jtree, cov(tgdata), n=nrow(tgdata), p.value, drop = TRUE)
  wmat <- skel$amat#gives adjacency matrix of the learned skeleton via LCD algorithm
  vset <- rownames(wmat)
  sep.pairs <- skel$sep.pairs
  cliques <- tg.jtree@cliques
  n.clique <- length(cliques)
  ###print(wmat)
  ###print(tg.jtree)
  #################################################################
  #################################################################
  #finding chordless cycles 
  ccycles<-ChordlessCycles(wmat)
  lcc<-length(ccycles)
  #finding chordless cycles of length 3
  ccyclesOf3<-list()
  vtemp<-c()
  if(lcc>0){
    for (i in 1:lcc){
      if(length(ccycles[[i]])==3){
        vtemp<-c(vtemp,i)
      }
    }
  }
  if(length(vtemp)>0){
    ccyclesOf3<-ccycles[vtemp]
  }
  lccOf3<-length(ccyclesOf3)
  ########################################################################
  ########## Finding possible parent-neighbors nodes for a given node in 
  ########## a chordless cycle ###########################################
  ########################################################################
  `findPNs`<-function(v,ccycle){
    vParentsNeighbors<-c()
    for (i in 1:length(ccycle)) {
      if(wmat[v,ccycle[i]]==5){
        vParentsNeighbors<-c(vParentsNeighbors,ccycle[i])
      }
    }
    return(vParentsNeighbors)
  }
  # ##########################################################################################################
  # ##########################################################################################################
  # ###A function that determines the type of chordless cycles: with or without blocks########################
  # ### if the cycle is chordless and without blocks, the function returns 1 o.w., returns 0 #################
  # ##########################################################################################################
  # ##########################################################################################################
  # ccycles_type<-function(vect,amat){
  #   cu<-cd<-co<-0
  #   l<-length(vect)
  #   for (i in 1:(l-1)) {
  #     if(amat[vect[i],vect[i+1]]==1 &&  amat[vect[i+1],vect[i]]==1){
  #       cu<-cu+1
  #     }
  #   }
  #   if(amat[vect[l],vect[1]]==1 &&  amat[vect[1],vect[l]]==1){
  #     cu<-cu+1
  #   }
  #   if(cu==l){
  #     #undirected cycle
  #     return(1)
  #   }else{
  #     return(0)
  #   }
  # }
  # ########################################################################################
  # ########## A function that returns the separation set of two given variables ###########
  # ########################################################################################
  `AsepCD`<-function(u,v,w){
    pair <- new("sep.pair", u=vset[u], v=vset[v])
    a <- lapply(sep.pairs, function(x) all.equal(x,pair))
    ###message(a," for ",u," and ", v)
    if(length((which(a==TRUE)))>0){
      sep<-sep.pairs[[which(a==TRUE)[1]]]@s
      if(is.element(vset[w], sep)){
        return(1)
      }else{
        return(0)
      }
    }else{
      return(2)
    }
  }
  
  ########################################################################################
  ########################################################################################
  ########## Then apply the rules: 1-4 ###################################################
  ########################################################################################
  change<-TRUE
  while(change){
    change<-FALSE
    ########################################################################
    ###### Rule 1 ###### A o--o B o--o C & B\not\in S_AC -> A +-o B o-+ C
    ########################################################################
    for(i in 1:n.clique){
      cvset <- cliques[[i]]$vset
      ##print(cvset)
      p <- length(cvset)
      if(p > 2){
        for(j in 1:(p-1)){
          for(k in (j+1):p){
            for(l in 1:p){
              u <- cvset[j]
              v <- cvset[k]
              w <- cvset[l]
              ####Then check A and C are not adjacent and B\not\in S_AC
              if((wmat[u,v]+wmat[v,u] == 0) &&
                 (wmat[w,u] != 5 || wmat[w,v] != 5) &&
                 (wmat[u,w]+wmat[w,u] > 0) && (wmat[w,v]+wmat[v,w] > 0)){
                pair <- new("sep.pair", u=u, v=v)
                a <- lapply(sep.pairs, function(x) all.equal(x,pair))
                ###message(a," for ",u," and ", v)
                sep <- sep.pairs[[which(a==TRUE)[1]]]@s
                ###message(sep," for ",u," and ", v)
                if(!is.element(w, sep)){
                  ###then block from B to A and C
                  wmat[w,u] <-5
                  wmat[w,v] <-5
                  change<-TRUE
                  ##message(w ," is a sudo-collider for ", u ," and ", v," rule 1")
                  ##print(wmat)
                  break
                }
              }#second if
              #count<- count + 1
            }#fourth for loop
          }#third for loop
        }#second for loop
      }#first if
      ###message("counter= ",count)
    }#first for loop
    #end of rule 1
    ########################################################################
    ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
    ########################################################################
    for (i in 1:length(V)) {
      ###message("change= ",change)
      bAdjs<-c()#parent, neighbor, or children of B
      bParentsOrNeighbors<-c()
      for (j in 1:length(V)) {
        if(wmat[j,i]+wmat[i,j]>0){
          bAdjs<-c(bAdjs,j)
        }
        if(wmat[i,j]==5){
          bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
        }
      }
      ###message(i," bParNeighbors: ",bParentsOrNeighbors)
      ####Then check A and C are not adjacent and B\in S_AC
      if(length(bAdjs)>1 && length(bParentsOrNeighbors)>0){
        for (k in 1:length(bParentsOrNeighbors)) {
          for(l in 1:length(bAdjs)){
            aIndex<-bParentsOrNeighbors[k]
            cIndex<-bAdjs[l]
            ##if A,C are not adjacents and B\in S_AC
            if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && (wmat[aIndex,cIndex]+ wmat[cIndex,aIndex]==0) && (AsepCD(aIndex,cIndex,i)==1 ||
                                                                                                            AsepCD(aIndex,cIndex,i)==2)){
              # ##message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
              # ##print(skel@sepset[[aIndex]][[cIndex]])
              # ##print(skel@sepset[[cIndex]][[aIndex]])
              ###then block from C to B
              wmat[cIndex,i]<-5
              ##message("b: ",vset[i]," "," c: ",vset[cIndex]," a: ",vset[aIndex]," rule 2")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 2
    ########################################################################
    ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
    ########################################################################
    if(lcc>0){
      for (i in 1:lcc) {
        cc<-ccycles[[i]]
        ###print(cc)
        #cc[j] plays the role of A
        for (j in 1:length(cc)) {
          aAdjs<-c()
          for (k in 1:length(cc)) {
            #aAdjs<-c()
            if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
              aAdjs<-c(aAdjs,cc[k])
            }
          }
          ###message(cc[j]," as A: ",aAdjs)
          if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
            #aAdjs[1] is a candidate for being B
            bPN<-findPNs(aAdjs[1],cc)
            if(length(bPN)==1){
              wmat[aAdjs[1],cc[j]]<-5
              ##message("aAdjs[1] ",aAdjs[1]," cc[j] ",cc[j], "rule 3a")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
          if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
            #aAdjs[2] is a candidate for being B
            bPN<-findPNs(aAdjs[2],cc)
            if(length(bPN)==1){
              wmat[aAdjs[2],cc[j]]<-5
              ##message("aAdjs[1] ",aAdjs[1]," cc[j] ",cc[j], "rule 3b")
              ##print(wmat)
              change<-TRUE
              break
            }
          }
        }
      }
    }#end of rule 3
    ###Now, since rule 1, 2, and 3 are faster than rule 4, apply these until no change, then go on to rule 4
    if (change){
      next
    }
    
    ########################################################################
    ###### Rule 4 ##########################################################
    ########################################################################
    
    if(lccOf3>1){
      for (i in 1:((lccOf3)-1)) {
        for (j in (i+1):lccOf3) {
          abIndex<-intersect(ccyclesOf3[[i]],ccyclesOf3[[j]])
          cIndex<-setdiff(ccyclesOf3[[i]],abIndex)
          dIndex<-setdiff(ccyclesOf3[[j]],abIndex)
          if((length(abIndex)==2) && (wmat[cIndex,dIndex]+wmat[dIndex,cIndex]==0)){
            for (turn in 1:2) {
              aIndex <- if(turn == 1) abIndex[1] else abIndex[2]
              bIndex <- if(turn == 1) abIndex[2] else abIndex[1]
              ##if C,D are not adjacents and A\in S_CD and if B to C,D is blocked and B to A is not blocked
              if(wmat[bIndex,aIndex]!=5 && wmat[aIndex,bIndex]!=0 && wmat[bIndex,cIndex]==5 && wmat[bIndex,dIndex]==5 &&
                 (AsepCD(cIndex,dIndex,aIndex)==1 || AsepCD(cIndex,dIndex,aIndex)==2)){
                #then block B to A i.e., A +- B
                wmat[bIndex,aIndex]<-5
                ##message("bIndex ",vset[bIndex]," aIndex ",vset[aIndex]," cIndex ",vset[cIndex]," dIndex ",vset[dIndex])
                ##print(wmat)
                change<-TRUE
                break
              }
            }
          }
        }
      }
    }##end of rule 4
    
    # for (i in 1:length(V)) {
    #   aAdjs<-c()
    #   for (j in 1:length(V)) {
    #     if(wmat[j,i]+wmat[i,j]>0){
    #       aAdjs<-c(aAdjs,j)
    #     }
    #   }
    #   if(length(aAdjs)>=3){
    #     for (k in 1:length(aAdjs)) {
    #       for (l in 1:length(aAdjs)) {
    #         cIndex<-aAdjs[k]
    #         dIndex<-aAdjs[l]
    #         ##if C,D are not adjacents and A\in S_CD
    #         if((cIndex!=dIndex) && !((wmat[i,cIndex]+wmat[cIndex,i]==2 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,dIndex]+wmat[dIndex,i]==2 && wmat[i,cIndex]==5 && wmat[cIndex,i]==1) ||
    #                                 (wmat[i,cIndex]==5 && wmat[cIndex,i]==1 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,cIndex]+wmat[cIndex,i]==10 && wmat[i,dIndex]==5 && wmat[dIndex,i]==1) ||
    #                                 (wmat[i,dIndex]+wmat[dIndex,i]==10 && wmat[i,cIndex]==5 && wmat[cIndex,i]==1))){
    #           for (r in 1:length(aAdjs)) {
    #             bIndex<-aAdjs[r]
    #             #and if B to C,D is blocked and B to A is not blocked
    #             if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
    #               #then block B to A i.e., A +- B
    #               wmat[bIndex,i]<-5
    #               ##message("bIndex ",vset[bIndex]," A ",vset[i]," cIndex ",vset[cIndex]," dIndex ",vset[dIndex])
    #               ##print(aAdjs)
    #               ##print(wmat)
    #               change<-TRUE
    #               #break
    #             }
    #           }
    #         }
    #       }
    #     }
    #   }
    # }##end of rule 4
  }
  # ##################################################################################################
  #################################################################################################### 
  # ################### This part of code can be used to obtain the Essential graph ##################
  # ################### of the Markov equivalence class of G and not just a CG in class ##############
  # ##### Learning marginal AMP chain graphs under faithfulness revisited (Pena & Gomez-Olmedo, 2016)#
  # ##################################################################################################
  # ##################################################################################################
  # ##################################################################################################
  # #  Replace every edge A - B in every cycle that is of length greater than three, chordless,####### 
  # #  and without blocks with A +-+ B ###############################################################
  # ##################################################################################################
  #   if(lcc>0){
  #     for (i in 1:lcc) {
  #       c<-ccycles[i]
  #       l<-length(c)
  #       if((length(c)>3) && (ccycles_type(c,wmat)==1)){
  #         for (j in 1:(l-1)) {
  #           wmat[c[j],c[j+1]]<-wmat[c[j+1],c[j]]<-5
  #         }
  #         wmat[c[1],c[l]]<-wmat[c[l],c[1]]<-5
  #       }
  #     }
  #   }
  # ##################################################################################################
  # ##################################################################################################
  # ##  Apply the rules R2-R4 while possible #########################################################
  # ##################################################################################################
  # ##################################################################################################
  #   change<-TRUE
  #   while(change){
  #     change<-FALSE
  #     ########################################################################
  #     ###### Rule 2 ####### A +-o B o--o C & B\in S_AC -> A +-o B +-o C
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       ###message("change= ",change)
  #       bAdjs<-c()#parent, neighbor, or children of B
  #       bParentsOrNeighbors<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           bAdjs<-c(bAdjs,j)
  #         }
  #         if(wmat[i,j]==5){
  #           bParentsOrNeighbors<-c(bParentsOrNeighbors,j)
  #         }
  #       }
  #       ###message(i," bParNeighbors: ",bParentsOrNeighbors)
  #       ####Then check A and C are not adjacent and B\in S_AC
  #       if(length(bAdjs)>0 && length(bParentsOrNeighbors)>0){
  #         for (k in 1:length(bParentsOrNeighbors)) {
  #           for(l in 1:length(bAdjs)){
  #             aIndex<-bParentsOrNeighbors[k]
  #             cIndex<-bAdjs[l]
  #             ##if A,C are not adjacents and B\in S_AC
  #             if((cIndex!=aIndex) && (wmat[cIndex,i]!=5) && wmat[aIndex,cIndex]==0 && wmat[cIndex,aIndex]==0 && is.element(vset[i],sepset(aIndex,cIndex))){
  #               # ##message("i= ",i," , aIndex= ",aIndex," , cIndex= ",cIndex)
  #               # ##print(skel@sepset[[aIndex]][[cIndex]])
  #               # ##print(skel@sepset[[cIndex]][[aIndex]])
  #               ###then block from C to B
  #               wmat[cIndex,i]<-5
  #               ###message("b: ",vset[i]," ","neighbor: ",vset[nbIndex])
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 2
  #     ########################################################################
  #     ###### Rule 3 ####### B o--o A +-o --- +-o B  -> B o-+ A +-o --- +-o B
  #     ########################################################################
  #     if(lcc>0){
  #       for (i in 1:lcc) {
  #         cc<-ccycles[[i]]
  #         ###print(cc)
  #         #cc[j] plays the role of A
  #         for (j in 1:length(cc)) {
  #           aAdjs<-c()
  #           for (k in 1:length(cc)) {
  #             #aAdjs<-c()
  #             if(wmat[cc[j],cc[k]]+wmat[cc[k],cc[j]]>0){
  #               aAdjs<-c(aAdjs,cc[k])
  #             }
  #           }
  #           ###message(cc[j]," as A: ",aAdjs)
  #           if((wmat[aAdjs[1],cc[j]]!=5) && (wmat[aAdjs[2],cc[j]]==5)){
  #             #aAdjs[1] is a candidate for being B
  #             bPN<-findPNs(aAdjs[1],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[1],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #           if((wmat[aAdjs[1],cc[j]]==5) && (wmat[aAdjs[2],cc[j]]!=5)){
  #             #aAdjs[2] is a candidate for being B
  #             bPN<-findPNs(aAdjs[2],cc)
  #             if(length(bPN)==1){
  #               wmat[aAdjs[2],cc[j]]<-5
  #               change<-TRUE
  #               break
  #             }
  #           }
  #         }
  #       }
  #     }#end of rule 3
  #     
  #     ########################################################################
  #     ###### Rule 4 ##########################################################
  #     ########################################################################
  #     for (i in 1:length(V)) {
  #       aAdjs<-c()
  #       for (j in 1:length(V)) {
  #         if(wmat[j,i]+wmat[i,j]>0){
  #           aAdjs<-c(aAdjs,j)
  #         }
  #       }
  #       if(length(aAdjs)>=3){
  #         for (k in 1:length(aAdjs)) {
  #           for (l in 1:length(aAdjs)) {
  #             cIndex<-aAdjs[k]
  #             dIndex<-aAdjs[l]
  #             ##if C,D are not adjacents and A\in S_CD
  #             if((cIndex!=dIndex) && is.element(vset[i],sepset(dIndex,cIndex))){
  #               for (r in 1:length(aAdjs)) {
  #                 bIndex<-aAdjs[r]
  #                 #and if B to C,D is blocked and B to A is not blocked
  #                 if((wmat[bIndex,i]!=5) && (wmat[bIndex,cIndex]==5) && (wmat[bIndex,dIndex]==5)){
  #                   #then block B to A i.e., A +- B
  #                   wmat[bIndex,i]<-5
  #                   change<-TRUE
  #                   break
  #                 }
  #               }
  #             }
  #           }
  #         }
  #       }
  #     }##end of rule 4
  #   } 
  ##################################################################################################
  ################## Replace every edge A +- B with A -> B #########################################
  ################## Replace every edge A +-+ B with A - B #########################################
  ##################################################################################################
  ###print(wmat)
  #wmat<-t(wmat)
  for (i in 1:(length(V)-1)) {
    for (j in (i+1):length(V)) {
      if(wmat[i,j]==5 && wmat[j,i]==5){
        wmat[i,j]<-wmat[j,i]<-1
      }
      if(wmat[i,j]==5 && wmat[j,i]==1){
        wmat[i,j]<-0
        #wmat[j,i]<-0
      }
      if(wmat[j,i]==5 && wmat[i,j]==1){
        wmat[j,i]<-0
        #wmat[i,j]<-0
      }
    }
  }
  return(wmat)
}
#################################################################################################################
#############################################################################################
######## Implementation of Algorithm 3 Pseudo-code for the DeflaggingProcedure (G) ##########
###################### (Roverato & Studeny, 2006) ###########################################
#############################################################################################
`Largest_DeflaggedAMPCG`<-function(amat){
  n<-nrow(amat)
  #############################################################################################
  ######## Implementation of Algorithm 1 Pseudo-code for the LabelingAlgorithm (G) ############
  ###################### (Roverato & Studeny, 2006) ###########################################
  `LabelingAlgorithm`<-function(amat){
    #replacing every line 'a-b' in G by 'a o-o b' i.e., line is free at a & b.
    # 4 is used for free side and 5 for blocked side #
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if(amat[i,j]==1 && amat[j,i]==1){
          amat[i,j]<-amat[j,i]<-4
        }
      }
    }
    #print(amat)
    possible<-TRUE
    while(possible){
      possible<-FALSE
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          for (k in 1:n) {
            if(i!=k && j!=k && amat[i,j]==0 && amat[j,i]==0 && amat[i,k]==1 && amat[k,i]==0 && (amat[j,k]==4 || amat[j,k]==5) && amat[k,j]==4){
              #message(i," a:",k," -x ",j)
              ### i -> k -o j is converted to i -> k -x j i.e., free side is blocked ###
              amat[k,j]<-5
              possible<-TRUE
              break
            }
            if(i!=k && j!=k && amat[i,j]==0 && amat[j,i]==0 && (amat[i,k]==4 || amat[i,k]==5) && amat[k,i]==4 && amat[j,k]==1 && amat[k,j]==0){
              #message(j," b:",k," -x ",i)
              ### j -> k -o i is converted to j -> k -x i i.e., free side is blocked ###
              amat[k,i]<-5
              possible<-TRUE
              break
            }
            if(i!=k && j!=k && amat[i,j]==0 && amat[j,i]==0 && (amat[i,k]+amat[k,i]>1) && (amat[j,k]==4) && (amat[k,j]==4 || amat[k,j]==5)){
              #message(i," c:",j," -x ",k)
              ### i - k o- j is converted to i - k x- j i.e., free side is blocked ###
              amat[j,k]<-5
              possible<-TRUE
              break
            }
            if(i!=k && j!=k && amat[i,j]==0 && amat[j,i]==0 && (amat[k,i]==4 || amat[k,i]==5) && amat[i,k]==4 && (amat[j,k]+amat[k,j]>1)){
              #message(j," d:",i," -x ",k)
              ### j - k o- i is converted to j - k x- i i.e., free side is blocked ###
              amat[i,k]<-5
              possible<-TRUE
              break
            }
            if(i!=k && j!=k && amat[i,j]==4 && (amat[j,i]==4 || amat[j,i]==5) && (amat[i,k]==5) && (amat[k,i]==4 || amat[k,i]==5) && (amat[j,k]==5 || amat[j,k]==4) && (amat[k,j]==5)){
              #message(k," e:",i," -x ",j)
              ### i -x k -x j o- i is converted to i -x k -x j x- i i.e., free side is blocked ###
              amat[i,j]<-5
              possible<-TRUE
              break
            }
            if(i!=k && j!=k && (amat[i,j]==4 || amat[i,j]==5) && amat[j,i]==4 && (amat[k,j]==4 || amat[k,j]==5) && amat[j,k]==5 && (amat[i,k]==4 || amat[i,k]==5) && (amat[k,i]==5)){
              #message(k," f:",j," -x ",i)
              ### j -x k -x i o- j is converted to j -x k -x i x- j i.e., free side is blocked ###
              amat[j,i]<-5
              possible<-TRUE
              break
            }
          }
        }
      }
    }
    return(amat)
  }
  #############################################################################################
  ######## Implementation of Algorithm 2 Pseudo-code for the DirectingAlgorithm (G^l)##########
  ###################### (Roverato & Studeny, 2006) ###########################################
  `DirectingAlgorithm`<-function(amat){
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if(amat[i,j]==4 && amat[j,i]==5){
          ## i x-o j is converted to i -> j
          amat[i,j]<-1
          amat[j,i]<-0
        }
        if(amat[i,j]==5 && amat[j,i]==4){
          ## j x-o i is converted to j -> i
          amat[j,i]<-1
          amat[i,j]<-0
        }
      }
    }
    ### unlabel the resulting graph
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        if((amat[i,j]==4 && amat[j,i]==4) || (amat[i,j]==5 && amat[j,i]==5) || (amat[i,j]==4 && amat[j,i]==5) || (amat[i,j]==5 && amat[j,i]==4)){
          amat[i,j]<-amat[j,i]<-1
        }
      }
    }
    return(amat)
  }
  ###############################################################################################
  ###############################################################################################
  temp<-amat
  temp[which(temp!=0)]<-0
  while(!identical(temp,amat)){
    temp<-amat
    amat<-LabelingAlgorithm(amat)
    amat<-DirectingAlgorithm(amat)
  }
  return(amat)
}
#####################################################################################################################################################