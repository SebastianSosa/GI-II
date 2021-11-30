library(ggplot2)
library(ggpubr)
library(igraph)
library(sna)
library(ANTs)
########################################################
# SRI for directed behaviours----------------------
directed.sri <- function(df, scan, actor = 'Actor', receiver = 'Receiver', weigth = "weigth",
                         method = 'sri', ynull = FALSE , num.ids = T){
  col_w <- ANTs:::df.col.findId(df, weigth)
  ### Get index of column with the scan
  col_scan <- ANTs:::df.col.findId(df, scan)
  if (length(col_scan) > 1) {
    df$scan <- apply(df[, col_scan ], 1, paste, collapse = "_")
  } else {
    df$scan <- df[, col_scan]
  }
  df$scan <- as.character(df$scan)
  
  ### Get all individuals' ids
  col_actor <- ANTs:::df.col.findId(df, actor)
  col_receiver <- ANTs:::df.col.findId(df, receiver)
  ids = unique(c(df[,col_actor], df[,col_receiver]))
  
  #Yab = matrix(0, ncol = length(ids), nrow = length(ids))
  #colnames(Yab) = rownames(Yab) = ids
  
  # Xab is the total number of interactions from a to b
  Xab = df.to.mat(df, col_actor, col_receiver, col_w, num.ids = T)
  Yab = Xab
  Yab[Yab > 0] = 0
  Ynull = Yab
  
  # Xba is the total number of interactions from b to a
  Xba = t(Xab)
  
  # Ya is the total number of interactions emmitted by a to other individuals (i.e. outstrength - x)
  Ya = matrix(met.outstrength(Xab), length(ids), length(ids), byrow = F)
  colnames(Ya) = rownames(Ya) = colnames(Xab)
  Ya = Ya - Xab
  Ya[Xab == 0] = 0
  
  # Yb is the total number of interactions received by b that are not from a (i.a instrength - x)
  Yb = matrix(met.instrength(Xab), length(ids), length(ids), byrow = T)
  colnames(Yb) = rownames(Yb) = colnames(Xab)
  Yb = Yb - Xab
  Yb[Xab == 0] = 0
  
  #Yab is the number of times a and b have been observed during the same scan but didn't interact.
  # Ynull is the number of times whereas neiher a and b have been observed.
  # To compute Yab and Ynull we will check on each scan whereas a & b have been observed together or not
  
  tmp = split(df, df[,col_scan])
  for (a in 1:length(tmp)) {
    # For Yab, creating a matrix with only individuals observed for the specific scan.
    # sym = TRUE as we want all individuals interactions in the scan
    mtmp = df.to.mat(tmp[[a]], col_actor, col_receiver, sym = T , num.ids = T) 
    #mtmp[upper.tri(mtmp, diag = TRUE)] = 0 # we remove lower triangle to avoid couting twice individual presence and as matrix fill through giving behaviours
    
    # For Ynull, Individuals prensent in mtmp are oberseved so we can extract non observed individuals
    ids.tmp = ids[!ids %in% colnames(mtmp)]
    Ynull[rownames(Ynull), colnames(Ynull) %in% ids.tmp] = Ynull[rownames(Ynull), colnames(Ynull) %in% ids.tmp] + 1
    
    if(nrow(mtmp) != 0){
      
      # converting interactions into NA
      mtmp[mtmp > 0 ] = NA
      # Converting non interactions in 1
      mtmp[!is.na(mtmp)] = 1
      # Diag = 0
      diag(mtmp) = 0
      # Converting interaction into 0
      mtmp[is.na(mtmp)] = 0
      
      # Matching colnames mtmp with matrix Yab and Ynull
      rowmatch <- match(rownames(mtmp), rownames(Yab))
      colmatch <- match(colnames(mtmp), colnames(Yab))
      
      Yab[rowmatch, colmatch] <- Yab[rowmatch, colmatch] + mtmp
    }
  }
  diag(Yab) = 0
  Yab = Yab[match(colnames(Xab), colnames(Yab)),match(colnames(Xab), colnames(Yab))]
  
  if(method == 'sri'){
    if(ynull){
      result = (Xab/(Xab + Xba + Ya + Yb + Yab + Ynull))
    }else{
      result = (Xab/(Xab + Xba + Ya + Yb + Yab))
    }
    
    result[is.na(result)] = 0
    return(result)
  }else{
    "Currently not available."
  }
  return(result)
}

### BASIC SIMULATION FUNCTIONS
############################################################################################################
######### Modification  1 : Global index approach
### function to create a network (matrix n x n) from data collection (obs) and focal id using simple ratio index (sri) 
make_network <- function(obs, focal.id) {
  # Creating data from simulation to allow ISI computation
  #obs = obs[-which(rowSums(obs) == 0),]
  df = NULL
  colnames(obs) = 1:ncol(obs)
  rownames(obs) = 1:nrow(obs)
  for(a in 1:nrow(obs)){
    ids.associated = names(which(obs[a, ] != 0))
    if(length(ids.associated) == 0){next()}
    if(focal.id[a] %in% ids.associated){stop()}
    d = data.frame("ego" = focal.id[a], "alter" = ids.associated, 'scan' = a)
    df = rbind(df, d)
  }
  df$weigth = 1
  m = directed.sri(df, scan = "scan", actor = "ego" , receiver = "alter", weigth = "weigth", ynull = FALSE)
  return(m)
}
### function to generate pre-network permutations (swaps of individuals between focals)
# No modifications
rand_network2 <- function(obs.p, focal.id, n.perm,n_focals) {
  N <- ncol(obs.p)
  networks_rand <- array(0, c(n.perm,N,N))
  for (i in 1:n.perm) {
    # first randomly select two focal observations
    repeat {
      o <- 1:n_focals
      a <- sample(o,1)
      b <- sample(o[-a],1)
      
      # check if these are different individuals and they have associates
      if ((focal.id[a] != focal.id[b]) & (sum(obs.p[a,])>0) & (sum(obs.p[b,])>0)) {
        # next select two associates to swap
        d <- sample(which(obs.p[a,] > 0),1)
        e <- sample(which(obs.p[b,] > 0),1)
        
        # check they do not occur in the other focal
        if ((obs.p[a,e] == 0) & obs.p[b,d] == 0) {
          
          # now check we have 4 distinct individuals, otherwise repeat this process
          if (!(d %in% c(focal.id[a], focal.id[b], e)) & !(e %in% c(focal.id[a], focal.id[b], d))) {
            break;
          }
        }
      }
    }
    
    # swap individuals
    obs.p[a,d] <- 0
    obs.p[b,d] <- 1
    obs.p[b,e] <- 0
    obs.p[a,e] <- 1
    # caculate network
    networks_rand[i,,] <- make_network(obs.p,focal.id)
  }
  return(networks_rand)
}

### Function to allocate number of observations to groups
# No modifications
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

### MAIN SIMULATION FUNCTION ###
### Arguments ###
# GS --> numeric argument indicating group size
# ObsBia --> numeric argument indicating the degree of observation bias [0.5-1.0]
# FemPhenotypeBias --> boolean argument indicating whether a phenotype bias is present among females
# nfocals --> numeric argument indicating number of focal samples
# N.perm --> numeric argument indicating number of permutations 
Simulation<-function(GS,ObsBias,FemSexRatio,FemPhenotypeBias,nfocals,N.Perm)
{
  # Set parameters
  N <- GS
  n_focals <- nfocals
  # Generate nodes
  NumFem<-round(GS * FemSexRatio)
  NumMal<-GS - NumFem
  Sex<-c(rep("F",NumFem),rep("M",NumMal))
  Sex<-sample(Sex,GS,replace=F)
  ids <- data.frame(ID=1:(N),SEX=Sex)
  # Generate a distribution of group sizes
  group_size <- sample(c(1:(N/2)),n_focals,replace=TRUE)
  # Create blank observation matrix
  obs <- matrix(0,nrow=n_focals,ncol=N)
  ## set number of observations of an individual in a group per individual
  ids$OBS <- rand_vect(N,sum(group_size),pos.only=TRUE)
  ## Variables to Allocate individuals to groups, 
  GroupID<-c(1:n_focals)
  group_size.tmp <- group_size
  # IF Fem phenotype is stronger than males, start with males so that they end up in smaller groups
  if(FemPhenotypeBias == T)
  {
    which.males <- which(ids$SEX=="M")
    which.females <- which(ids$SEX=="F")
    for (i in which.males) 
    {
      g <- sample(GroupID[which(group_size.tmp>0)],ids$OBS[i])
      group_size.tmp[g] <- group_size.tmp[g]-1
      obs[g,i] <- 1
    }
    for (i in which.females) 
    {
      if ((sum(group_size.tmp>0) < ids$OBS[i])) 
      {
        Needed<-ids$OBS[i]-(sum(group_size.tmp>0))
        group.tmp<-group_size
        group.tmp[group_size.tmp>0]=0
        BiggestGroups<-sort(group.tmp,decreasing = T,index.return=T)$ix
        ExtraGroups<-BiggestGroups[1:Needed]
        g<-c(GroupID[which(group_size.tmp>0)],ExtraGroups)
      }else 
      {
        g<-sample(GroupID[which(group_size.tmp>0)],ids$OBS[i])
      }
      group_size.tmp[g] <- group_size.tmp[g]-1
      obs[g,i] <- 1
    }
  }else # IF Fem phenotype is equal to males, allocate indivdiuals to groups at random
  {
    Inds<-c(1:GS)
    for (. in 1:GS) 
    {
      id<-Inds[1]
      if(length(Inds)>1){id<-sample(Inds,1)}
      Inds<-Inds[-which(Inds==id)]
      if ((sum(group_size.tmp>0) < ids$OBS[id])) 
      {
        Needed<-ids$OBS[id]-(sum(group_size.tmp>0))
        Fullgroups<-which(group_size.tmp==0)
        ExtraGroups<-sample(Fullgroups,Needed,replace=F) 
        g<-c(GroupID[which(group_size.tmp>0)],ExtraGroups)
      }else 
      {
        g<-sample(GroupID[which(group_size.tmp>0)],ids$OBS[id])
      }
      group_size.tmp[g] <- group_size.tmp[g]-1
      obs[g,id] <- 1
    }
    
  }
  # Select a focal individual from each group
  focal.id <- apply(obs,1,function(x) { sample(which(x==1),1)})
  
  # Now remove cases where individuals occur in a group for which they are focal
  obs[cbind(1:n_focals,focal.id)] <- 0
  
  ## NOW DO NETWORK ANALYSIS ON THESE DATA
  # Calculate network
  Net.Ori <- make_network(obs,focal.id)
  sampling.effort = Net.Ori
  sampling.effort[sampling.effort>0] = 0
  for(a in 1:length(focal.id)){
    sampling.effort[focal.id[a],] = sampling.effort[focal.id[a],] + 1
    sampling.effort[,focal.id[a]] = sampling.effort[,focal.id[a]] + 1
  }
  diag(sampling.effort) = 0

  Net.Ori.corrected <- Net.Ori/sampling.effort
  diag(Net.Ori.corrected) = 0
  Net.Ori.corrected[is.nan(Net.Ori.corrected)] = 0
  Net.Ori.corrected[is.na(Net.Ori.corrected)] = 0
  Net.Ori.corrected[is.infinite(Net.Ori.corrected)] = 0
  # Remove some observations according to the degre of observation bias ObsBias
  # Generate probability of being observed (males=1,females=ObsBias)
  ids$OBS_PROB <- ObsBias
  ids$OBS_PROB[which(ids$SEX=="M")] <- 1
  
  # Remove observations from GBI
  obs.Bias <- obs
  for (i in 1:N) {
    obs.Bias[which(obs.Bias[,i] > 0),i] <- sample(c(0,1),sum(obs.Bias[,i]),replace=TRUE,prob=c(1-ids$OBS_PROB[i],ids$OBS_PROB[i]))
  }
  # Calculate new network
  Net.Bias <- make_network(obs.Bias,focal.id)
  sampling.effort = Net.Bias
  sampling.effort[sampling.effort>0] = 0
  for(a in 1:length(focal.id)){
    sampling.effort[focal.id[a],] = sampling.effort[focal.id[a],] + 1
    sampling.effort[,focal.id[a]] = sampling.effort[,focal.id[a]] + 1
  }
  diag(sampling.effort) = 0
  
  Net.Biais.corrected <- Net.Ori/sampling.effort
  diag(Net.Biais.corrected) = 0

  Net.Biais.corrected[is.nan(Net.Biais.corrected)] = 0
  Net.Biais.corrected[is.na(Net.Biais.corrected)] = 0
  Net.Biais.corrected[is.infinite(Net.Biais.corrected)] = 0
  # Calculate Strength
  ids$DEGREE <- rowSums(Net.Ori)
  ids$DEGREE.Corrected <- rowSums(Net.Ori)
  ids$DEGREE.Bias <- rowSums(Net.Bias)
  ids$DEGREE.Bias.Corrected <- rowSums(Net.Biais.corrected)
  print(summary(lm(DEGREE.Bias~SEX,data=ids)))
  
  obs.per.ind.Bias =   rep(0, nrow(ids))
  for (x in 1:nrow(ids)) {
    obs.per.ind.Bias[x] = length(which(focal.id %in% x))
  }
  ids$obs.bias = obs.per.ind.Bias
  
  p1 = ggplot(ids, aes(x = SEX, y = obs.bias, color = SEX))+geom_text(label = ids$ID)+ggtitle("Bias of observaiton")
  p2 = ggplot(ids, aes(x = SEX, y = DEGREE, color = SEX))+geom_text(label = ids$ID)+ggtitle("True relationship between strength and sex")
  p3 = ggplot(ids, aes(x = SEX, y = DEGREE.Corrected, color = SEX))+geom_text(label = ids$ID)+ggtitle("Biased relationship between strength and sex")
  p4 = ggplot(ids, aes(x = SEX, y = DEGREE.Bias.Corrected, color = SEX))+geom_text(label = ids$ID)+ggtitle("Corrected relationship between strength and sex")
  p5 = ggplot(ids, aes(x = DEGREE, y = DEGREE.Bias, color = SEX))+geom_text(label = ids$ID)+ggtitle("Correlation between true strength and biased")
  p6 = ggplot(ids, aes(x = DEGREE.Corrected, y = DEGREE.Bias.Corrected, color = SEX))+geom_text(label = ids$ID)+ggtitle("Correlation between true strength and corrected")
  
  print(ggarrange(p1, p2, p3, p4, p5, p6,  ncol = 3, nrow = 2, common.legend = T))
  
  ############################################################################################################
  ######### Modification  2 : Compute degree and eigenvector
  ids$alters <- met.outdegree(Net.Ori)
  alters.bias <- met.outdegree(Net.Biais.corrected)
  ids$alters.Bias <- (alters.bias)/ obs.per.ind.Bias
  if(any(is.infinite(ids$alters.Bias))){ids$alters.Bias[which(is.infinite(ids$alters.Bias))] = NA}
  
  
  ids$DEGREE.Bias =  ((ids$DEGREE.Bias.Corrected ))
  print(summary(lm(DEGREE.Bias~SEX,data=ids)))
  
  
  ids$eigen <- met.eigen(Net.Ori)
  ids$eigen.Bias <- ((met.eigen(Net.Biais.corrected)))
  
  ############################################################################################################
  ######### Modification 3 : Running simulations for degree, eigenvector to
  
  # Calculate effects
  coef.Ori <- coefficients(lm(DEGREE~SEX,data=ids))[2]
  coef.Bias <- coefficients(lm(DEGREE.Bias~SEX,data=ids))[2]
  cat("Bias coefficient: ", coef.Bias, "\n")
  if(coef.Bias > 0){warning("relationship is inverted")}
  cat("Amount of bias: ", ObsBias, "\n")
  coef.eigen.Bias <- coefficients(lm(eigen.Bias~SEX,data=ids))[2]
  coef.alters.Bias <- coefficients(lm(alters.Bias~SEX,data=ids))[2]
  ### Data permutations
  n.perm <- N.Perm
  coefs.Perm_Nodes = coefs.eigen.Perm_Nodes = coefs.alters.Perm_Nodes = NULL
  for(d in 1:N.Perm){
    coefs.Perm_Nodes[d] = summary(lm(data = ids, formula = DEGREE.Bias ~ sample(SEX)))$coefficients[2,1]
    coefs.eigen.Perm_Nodes[d] = summary(lm(data = ids, formula = eigen.Bias ~ sample(SEX)))$coefficients[2,1]
    coefs.alters.Perm_Nodes[d] = summary(lm(data = ids, formula = alters.Bias ~ sample(SEX)))$coefficients[2,1]
  }
  #if(sum(coef.Bias>coefs.Perm_Nodes) / n.perm > 0.05){stop()}
  ############################################################################################################
  ######### Modification  6 : Returning only p-values
  Result <- data.frame(#"Strength pre-network" = sum(coef.Bias>coefs_Perm) / n.perm, 
                       "Strength network" = sum(coef.Bias>coefs.Perm_Nodes) / n.perm,
                       "Strength parametric" = summary(lm(DEGREE.Bias~SEX,data=ids))$coefficients[2,4],
                       #"Eigen pre-network" = sum(coef.eigen.Bias>coefs_eigen_Perm) / n.perm, 
                       "Eigen network" = sum(coef.eigen.Bias>coefs.eigen.Perm_Nodes) / n.perm,
                       "Eigen parametric" = summary(lm(eigen.Bias~SEX,data=ids))$coefficients[2,4],
                       #"Alters pre-network" = sum(coef.alters.Bias>coefs_alters_Perm) / n.perm,
                       "Alters network" = sum(coef.alters.Bias>coefs.alters.Perm_Nodes) / n.perm,
                       "Alters parametric" = summary(lm(alters.Bias~SEX,data=ids))$coefficients[2,4])
  
  Result
}

##################################
# Latin hypercube sampling
###################
## Simulations with biases of observation
library(lhs)
NumCombinations<-500
VariablesToSample<-4
VarNames<-c("GroupSize",    ## Range 30-100
            "FEM.REMOVAL",  ## Range 0.5-1
            "FemSexRatio",  ## Range 0.2-0.8
            "FOCALS.NUM")   ## Range 100 - 2000
LHS<-randomLHS(NumCombinations,VariablesToSample)
Mat<-matrix(NA,nrow=NumCombinations,ncol=4)
Mat[,1]<-round((30 + (LHS[,1]*(100-10))),0)
Mat[,2]<-round(0.5 + (LHS[,2]*(1-0.5)),2)
Mat[,3]<-round(0.2 + (LHS[,3]*(0.8-0.2)),2)
Mat[,4]<-round((100 + (LHS[,4]*(2000-100))),0)
FemPhenotypeBias<-c(TRUE,FALSE)
nSim = 1
R =  NULL
for (a in 1:length(FemPhenotypeBias))
{
  for(b in 1:nrow(Mat))
  {
    Result<-list()
    result = NULL
    for(c in 1:nSim)
    {
      cat("Simulation: ", b, "\n")
      
      
      df = Simulation(
        GS = Mat[b,1],ObsBias = Mat[b,2], FemSexRatio = Mat[b,3],FemPhenotypeBias = FemPhenotypeBias[a], nfocals = Mat[b,4],
        N.Perm = 1000)
      
      df$GS = Mat[b,1]
      df$ObsBias = Mat[b,2]
      df$FemSexRatio = Mat[b,3]
      df$FemPhenotypeBias = FemPhenotypeBias[a]
      df$nfocals =Mat[b,4] 
      df
      R = rbind(R, df)
      
      cat("Parametric true positive rates for strength: ", sum(R[R$FemPhenotypeBias == T,]$Strength.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for strength: ", sum(R[R$FemPhenotypeBias == T,]$Strength.network <0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for strength: ", sum(R[R$FemPhenotypeBias == T,]$Strength.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      
      cat("Parametric true positive rates for eigenvector: ", sum(R[R$FemPhenotypeBias == T,]$Eigen.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for eigenvector: ", sum(R[R$FemPhenotypeBias == T,]$Eigen.network <0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for eigenvector: ", sum(R[R$FemPhenotypeBias == T,]$Eigen.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      
      cat("Parametric true positive rates for alters: ", sum(R[R$FemPhenotypeBias == T,]$Alters.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for alters: ", sum(R[R$FemPhenotypeBias == T,]$Alters.network <0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for alters: ", sum(R[R$FemPhenotypeBias == T,]$Alters.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == T,]), "\n")
      
      
      cat("Parametric false positive rates for strength: ", sum(R[R$FemPhenotypeBias == F,]$Strength.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation false positive rates for strength: ", sum(R[R$FemPhenotypeBias == F,]$Strength.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation false positive rates for strength: ", sum(R[R$FemPhenotypeBias == F,]$Strength.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      
      cat("Parametric false positive rates for eigenvector: ", sum(R[R$FemPhenotypeBias == F,]$Eigen.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation false positive rates for eigenvector: ", sum(R[R$FemPhenotypeBias == F,]$Eigen.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation false positive rates for eigenvector: ", sum(R[R$FemPhenotypeBias == F,]$Eigen.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      
      cat("Parametric false positive rates for alters: ", sum(R[R$FemPhenotypeBias == F,]$Alters.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation false positive rates for alters: ", sum(R[R$FemPhenotypeBias == F,]$Alters.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation false positive rates for alters: ", sum(R[R$FemPhenotypeBias == F,]$Alters.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]), "\n")
    }
  }
}


## Simulations without biases of observation
Mat[,2]<-rep(1,NumCombinations) # no obs bias keep constant to 1
R2 =  NULL
for (a in 1:length(FemPhenotypeBias))
{
  for(b in 1:nrow(Mat))
  {
    Result<-list()
    result = NULL
    for(c in 1:nSim)
    {
      cat(b, "\n")
      
      
      df = Simulation(
        GS = Mat[b,1],ObsBias = Mat[b,2], FemSexRatio = Mat[b,3],FemPhenotypeBias = FemPhenotypeBias[a], nfocals = Mat[b,4],
        N.Perm = 100)
      
      df$GS = Mat[b,1]
      df$ObsBias = Mat[b,2]
      df$FemSexRatio = Mat[b,3]
      df$FemPhenotypeBias = FemPhenotypeBias[a]
      df$nfocals =Mat[b,4] 
      R2 = rbind(R2, df)
      
      cat("Parametric true positive rates for strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      
      cat("Parametric true positive rates for eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      
      cat("Parametric true positive rates for alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      
      
      cat("Parametric false positive rates for strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation false positive rates for strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation false positive rates for strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      
      cat("Parametric false positive rates for eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation false positive rates for eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation false positive rates for eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      
      cat("Parametric false positive rates for alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation false positive rates for alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation false positive rates for alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
    }
  }
}
