library(ANTs)
library(ggplot2)
library(ggpubr)
### function to create a network (matrix n x n) from data collection (obs) and focal id using simple ratio index (sri)
############################################################################################################
######### Modification  1 (conceptual): Global index approach
make_network <- function(obs, focal.id) {
  N <- ncol(obs)
  network <- matrix(0,nrow=N, ncol=N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      xab <- sum(obs[which(focal.id %in% c(i,j)),c(i,j)]>0)
      ############################################################################################################
      ######### Modification  1.1 : Computing SRI where xab is the association between a and b,
      #########                     ya is the number of observations of a without b, and
      #########                     yb is the number of observations of b without a
      fa <- sum(focal.id == i)
      fb <- sum(focal.id == j)
      ya <- fa - xab
      yb <- fb - xab
      if(yb < 0){yb = 0}
      if(ya < 0){ya = 0}
      sri <- ((xab) /(xab + ya + yb ))

      if(!is.nan(sri)){
        network[i,j] <- sri
        network[j,i] <- sri
      }else{
        network[i,j] <- 0
        network[j,i] <- 0
      }
    }
  }
  return(network)
}

make_network.corrected <- function(obs, focal.id) {
  N <- ncol(obs)
  network <- matrix(0,nrow=N, ncol=N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      xab <- sum(obs[which(focal.id %in% c(i,j)),c(i,j)]>0)
      ############################################################################################################
      ######### Modification  1.2 : Computing SRI where xab is the association between a and b,
      #########                     ya is the number of observations of a without b, and
      #########                     yb is the number of observations of b without a
      fa <- sum(focal.id == i)
      fb <- sum(focal.id == j)
      ya <- fa - xab
      yb <- fb - xab
      if(yb < 0){yb = 0}
      if(ya < 0){ya = 0}
      sri <- ((xab) /(xab + ya + yb ))
      gi = fa + fb
      if(!is.nan(sri)){
        network[i,j] <- sri
        network[j,i] <- sri

      }else{
        network[i,j] <- 0
        network[j,i] <- 0
      }
    }
  }

  # Sampling effort for each individuals--------------------
  obs.per.ind.Bias =   rep(0, nrow(network))
  for (x in 1:ncol(network)) {
    obs.per.ind.Bias[x] = length(which(focal.id %in% x))
  }
  mean.obs = mean(obs.per.ind.Bias, na.rm = T)
  dif.mean = mean.obs - obs.per.ind.Bias
  dif.mean.obs = abs(dif.mean)
  #dif.mean = abs(dif.mean.obs)

  int = colSums(obs)
  mean.obs = mean(int, na.rm = T)
  dif.mean = mean.obs - int
  dif.mean = abs(dif.mean)

  # Select your diff.mean
  #if.mean = dif.mean/dif.mean.obs
  #dif.mean = dif.mean.obs # similar false negatives rates
  #dif.mean = int/obs.per.ind.Bias # 100% false negatives
  #dif.mean = int
  #dif.mean = obs.per.ind.Bias
  #dif.mean
  fa = matrix(0, ncol = nrow(network), nrow = nrow(network))
  for (a in 1:nrow(network)) {
    fa[a,] = dif.mean[a]
  }

  fb = t(fa)

  sampling.effort = log(fa + fb)
  diag(sampling.effort) = 0
  #range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  #sampling.effort = range01(sampling.effort)
  plot(as.vector(sampling.effort), as.vector(network))


  GI = (network/(sampling.effort))
  return(GI)
}

############################################################################################################

### function to generate pre-network permutations (swaps of individuals between focals)
# No modifications
rand_network <- function(obs.p, focal.id, n.perm,n_focals) {
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
# Data stream with GI
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
    networks_rand[i,,] <- make_network.corrected(obs.p,focal.id)
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

### MAIN SIMULATION FUNCTION#########
#' @param  GS numeric argument indicating group size
#' @param  ObsBia numeric argument indicating the degree of observation bias [0.5-1.0]
#' @param  FemPhenotypeBias boolean argument indicating whether a phenotype bias is present among females
#' @param  nfocals numeric argument indicating number of focal samples
#' @param  N.perm numeric argument indicating number of permutations
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
  Net.Ori.corrected <- make_network.corrected(obs,focal.id)
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
  Net.Biais <- make_network(obs.Bias,focal.id)
  Net.Biais.corrected <- make_network.corrected(obs.Bias,focal.id)
  Net.Biais.corrected[is.infinite(Net.Biais.corrected)] = 0
  Net.Biais.corrected[is.nan(Net.Biais.corrected)] = 0

  # Calculate Strength
  ids$DEGREE <- rowSums(Net.Ori)
  ids$DEGREE.Corrected <- rowSums(Net.Ori.corrected)

  ids$DEGREE.Bias <- rowSums(Net.Biais)
  ids$DEGREE.Bias.Corrected <- rowSums(Net.Biais.corrected)


  # Sampling effort for each individuals--------------------
  obs.per.ind.Bias =   rep(0, nrow(ids))
  for (x in 1:nrow(ids)) {
    obs.per.ind.Bias[x] = length(which(focal.id %in% x))
  }
  ids$obs.bias = obs.per.ind.Bias

  ############################################################################################################
  ######### Modification  2 (extension): Compute degree and eigenvector
  ids$alters <- met.degree(Net.Ori)
  ids$alters.Bias <- met.degree(Net.Biais)
  ids$alters.Bias.Corrected <- (ids$alters.Bias)/ obs.per.ind.Bias
  if(any(is.infinite(ids$alters.Bias.Corrected))){ids$alters.Bias.Corrected[which(is.infinite(ids$alters.Bias.Corrected))] = NA}

  ids$eigen <- met.eigen(Net.Ori)
  ids$eigen.Bias <- ((met.eigen(Net.Biais)))
  ids$eigen.Bias.Corrected <- ((met.eigen(Net.Biais.corrected)))

  ############################################################################################################
  ######### results visualization
  ############################################################################################################
  p1 = ggplot(ids, aes(x = SEX, y = DEGREE, group = SEX))+ geom_boxplot()+geom_point()
  p2 = ggplot(ids, aes(x = SEX, y = DEGREE.Bias, group = SEX))+ geom_boxplot()+geom_point()
  p3 = ggplot(ids, aes(x = SEX, y = DEGREE.Bias.Corrected, group = SEX))+ geom_boxplot()+geom_point()
  p4 = ggplot(ids, aes(x = SEX, y = DEGREE.Bias.Corrected2, group = SEX))+ geom_boxplot()+geom_point()

  p4 = ggplot(ids, aes(x = SEX, y = alters, group = SEX))+ geom_boxplot()+geom_point()
  p5 = ggplot(ids, aes(x = SEX, y = alters.Bias, group = SEX))+ geom_boxplot()+geom_point()
  p6 = ggplot(ids, aes(x = SEX, y = alters.Bias.Corrected, group = SEX))+ geom_boxplot()+geom_point()

  p7 = ggplot(ids, aes(x = SEX, y = eigen, group = SEX))+ geom_boxplot()+geom_point()
  p8 = ggplot(ids, aes(x = SEX, y = eigen.Bias, group = SEX))+ geom_boxplot()+geom_point()
  p9 = ggplot(ids, aes(x = SEX, y = eigen.Bias.Corrected, group = SEX))+ geom_boxplot()+geom_point()

  #p10 = ggplot(data=NULL, aes(x = ids$SEX, y= obs.per.ind.Bias, group = ids$SEX))+geom_boxplot()+geom_point()
  #print(ggarrange(p1, p2,p3, p4, p10, ncol = 3, nrow = 2))

  print(ggarrange(p1, p2,p3,p4,p5, p6, p7, p8, p9, ncol = 3, nrow = 3))


  ###### Data permutations ###############
  ############################################################################################################
  ######### Modification  5 : Node label with and without GI & double permutations with and without GI

  coefs.Perm_Nodes = coefs.eigen.Perm_Nodes = coefs.alters.Perm_Nodes = NULL
  coefs.Perm_Nodes_Corrected  = coefs.eigen.Perm_Nodes_Corrected  = coefs.alters.Perm_Nodes_Corrected  = NULL

  coefs.Perm_Nodes_dp = coefs.eigen.Perm_Nodes_dp = coefs.alters.Perm_Nodes_dp = NULL
  coefs.Perm_Nodes_dp_Corrected  = coefs.eigen.Perm_Nodes_dp_Corrected  = coefs.alters.Perm_Nodes_dp_Corrected  = NULL

  for(d in 1:N.Perm){
    # Node label
    coefs.Perm_Nodes[d] = summary(lm(data = ids, formula = DEGREE.Bias ~ sample(SEX)))$coefficients[2,1]
    coefs.eigen.Perm_Nodes[d] = summary(lm(data = ids, formula = eigen.Bias ~ sample(SEX)))$coefficients[2,1]
    coefs.alters.Perm_Nodes[d] = summary(lm(data = ids, formula = alters.Bias ~ sample(SEX)))$coefficients[2,1]

    coefs.Perm_Nodes_Corrected[d] = summary(lm(data = ids, formula = DEGREE.Bias.Corrected ~ sample(SEX)))$coefficients[2,1]
    coefs.eigen.Perm_Nodes_Corrected[d] = summary(lm(data = ids, formula = eigen.Bias.Corrected ~ sample(SEX)))$coefficients[2,1]
    coefs.alters.Perm_Nodes_Corrected[d] = summary(lm(data = ids, formula = alters.Bias.Corrected ~ sample(SEX)))$coefficients[2,1]
  }

  coef.Ori <- coefficients(lm(DEGREE~SEX,data=ids))[2]
  coef.Bias <- coefficients(lm(DEGREE.Bias~SEX,data=ids))[2]
  coef.Bias.Corrected  <- coefficients(lm(DEGREE.Bias.Corrected~SEX,data=ids))[2]

  print(summary(lm(DEGREE.Bias.Corrected~SEX,data=ids)))
  cat("Amount of bias: ", ObsBias, "\n")

  coef.eigen.Bias <- coefficients(lm(eigen.Bias~SEX,data=ids))[2]
  coef.eigen.Bias.Corrected <- coefficients(lm(eigen.Bias.Corrected~SEX,data=ids))[2]

  coef.alters.Bias  <- coefficients(lm(alters.Bias~SEX,data=ids))[2]
  coef.alters.Bias.Corrected  <- coefficients(lm(alters.Bias.Corrected~SEX,data=ids))[2]
  ############################################################################################################
  ######### Modification  6 : p-values permuted
  p1 = ANTs:::stat.p(c(coef.Bias, coefs.Perm_Nodes))
  p2 = ANTs:::stat.p(c(coef.eigen.Bias, coefs.eigen.Perm_Nodes))
  p3 = ANTs:::stat.p(c(coef.alters.Bias, coefs.alters.Perm_Nodes))

  p4 = ANTs:::stat.p(c(coef.Bias.Corrected, coefs.Perm_Nodes_Corrected))
  p5 = ANTs:::stat.p(c(coef.eigen.Bias.Corrected, coefs.eigen.Perm_Nodes_Corrected))
  p6 = ANTs:::stat.p(c(coef.alters.Bias.Corrected, coefs.alters.Perm_Nodes_Corrected))

  cat(p4,p5,p6, "\n")
  ############################################################################################################
  ######### Modification  6 : One-tailed parametric test
  s.degree = summary(lm(DEGREE.Bias~SEX,data=ids))
  s.eigen = summary(lm(eigen.Bias~SEX,data=ids))
  s.alters = summary(lm(alters.Bias~SEX,data=ids))

  p7 = pt(coef(s.degree)[,3], s.degree$df[2], lower = T)[2]
  p8 = pt(coef(s.eigen)[,3], s.eigen$df[2], lower = T)[2]
  p9 = pt(coef(s.alters)[,3], s.alters$df[2], lower = T)[2]

  p10 = pt(coef(s.degree)[,3], s.degree$df[2], lower = F)[2]
  p11 = pt(coef(s.eigen)[,3], s.eigen$df[2], lower = F)[2]
  p12 = pt(coef(s.alters)[,3], s.alters$df[2], lower = F)[2]

  p13 = s.degree$coefficients[2,4]
  p14 = s.eigen$coefficients[2,4]
  p15 = s.alters$coefficients[2,4]


  s.degree.Corrected  = summary(lm(DEGREE.Bias.Corrected~SEX,data=ids))
  s.eigen.Corrected  = summary(lm(eigen.Bias.Corrected~SEX,data=ids))
  s.alters.Corrected  = summary(lm(alters.Bias.Corrected~SEX,data=ids))

  p16 = pt(coef(s.degree.Corrected)[,3], s.degree.Corrected$df[2], lower = T)[2]
  p17 = pt(coef(s.eigen.Corrected)[,3], s.eigen.Corrected$df[2], lower = T)[2]
  p18 = pt(coef(s.alters.Corrected)[,3], s.alters.Corrected$df[2], lower = T)[2]

  p19 = pt(coef(s.degree.Corrected)[,3], s.degree.Corrected$df[2], lower = F)[2]
  p20 = pt(coef(s.eigen.Corrected)[,3], s.eigen.Corrected$df[2], lower = F)[2]
  p21 = pt(coef(s.alters.Corrected)[,3], s.alters.Corrected$df[2], lower = F)[2]

  p22 = s.degree.Corrected$coefficients[2,4]
  p23 = s.eigen.Corrected$coefficients[2,4]
  p24 = s.alters.Corrected$coefficients[2,4]

  ############################################################################################################
  ######### Modification  7 : Returning only p-values

  Result <- data.frame("Strength network" = p1,
                       "Strength parametric" = c(p7, p10,p13),


                       "Eigen network" = p2,
                       "Eigen parametric" =  c(p8, p11,p14),

                       "Alters network" = p3,
                       "Alters parametric" = c(p9, p12,p15),

                       "Strength network corrected" = p4,
                       "Strength parametric corrected" = c(p16, p19, p22),

                       "Eigen network corrected" = p5,
                       "Eigen parametric corrected" = c(p17, p20, p23),

                       "Alters network corrected" = p6,
                       "Alters parametric corrected" = c(p18, p21, p24)
                       )

  Result = t(Result)
  return(list(Result, coef.Bias.Corrected))
}

##################################
# Latin hypercube sampling
###################
## Simulations with biases of observation-------------------
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
nSim = a = b = c = 1
R =  NULL
for (a in a:length(FemPhenotypeBias))
{
  for(b in 1:nrow(Mat))
  {
    Result<-list()
    result = NULL
    for(c in 1:nSim)
    {
      cat("#################################################################################", '\n')
      cat("Simulation: ", b, "\n")


        d = Simulation(
          GS = Mat[b,1],ObsBias = Mat[b,2], FemSexRatio = Mat[b,3],FemPhenotypeBias = FemPhenotypeBias[a], nfocals = Mat[b,4],
          N.Perm = 1000)
        df = as.data.frame(d[[1]])
        df$GS = Mat[b,1]
        df$ObsBias = Mat[b,2]
        df$FemSexRatio = Mat[b,3]
        df$FemPhenotypeBias = FemPhenotypeBias[a]
        df$nfocals =Mat[b,4]
        df$test = rownames(df)
        df$sim = b
        df$coef = d[[2]]

      R = rbind(R, df)

      d1 = data.frame(
        "Biases" = rep(TRUE, 2),
        "Error Type" = rep("False negatives rates", 2),
        "GI" = rep(FALSE, 2),
        "two.side" = rep(TRUE, 2),
        "approches" = c("Parametric", "Nertwork permutations"),
        "strength" = c(sum(R[R$FemPhenotypeBias == T & R$test %in% "Strength.parametric",]$p_value_two_side>0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Strength.parametric",]),
                       sum(R[R$FemPhenotypeBias == T & R$test %in% "Strength.network",]$p_value_two_side >0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Strength.network",])),
        "eigenvector" = c(sum(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.parametric",]$p_value_two_side>0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.parametric",]),
                          sum(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.network",]$p_value_two_side >0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.network",])),
        "Alters" = c(sum(R[R$FemPhenotypeBias == T & R$test %in% "Alters.parametric",]$p_value_two_side>0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Alters.parametric",]),
                     sum(R[R$FemPhenotypeBias == T & R$test %in% "Alters.network",]$p_value_two_side >0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Alters.network",]))

      )


      d2 =  data.frame(
        "Biases" = rep(TRUE, 2),
        "Error Type" = rep("False negatives rates", 2),
        "GI" = rep(TRUE, 2),
        "two.side" = rep(TRUE, 2),
        "approches" = c("Parametric", "Nertwork permutations"),
        "strength" = c(sum(R[R$FemPhenotypeBias == T & R$test %in% "Strength.parametric.corrected",]$p_value_two_side>0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Strength.parametric.corrected",]),
                       sum(R[R$FemPhenotypeBias == T & R$test %in% "Strength.network.corrected",]$p_value_two_side >0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Strength.network.corrected",])),
        "eigenvector" = c(sum(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.parametric.corrected",]$p_value_two_side>0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.parametric.corrected",]),
                          sum(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.network.corrected",]$p_value_two_side >0.05)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Eigen.network.corrected",])),
        "Alters" = c(sum(R[R$FemPhenotypeBias == T & R$test %in% "Alters.parametric.corrected",]$p_value_two_side>0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Alters.parametric.corrected",]),
                     sum(R[R$FemPhenotypeBias == T & R$test %in% "Alters.network.corrected",]$p_value_two_side >0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T & R$test %in% "Alters.network.corrected",]))

      )

      d3 = data.frame(
        "Error Type" = rep("False positives rates", 2),
        "Biases" = rep(TRUE, 2),
        "GI" = rep(FALSE, 2),
        "two.side" = rep(TRUE, 2),
        "approches" = c("Parametric", "Nertwork permutations"),
        "strength" = c(sum(R[R$FemPhenotypeBias == F & R$test %in% "Strength.parametric",]$p_value_two_side<0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Strength.parametric",]),
                       sum(R[R$FemPhenotypeBias == F & R$test %in% "Strength.network",]$p_value_two_side <0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Strength.network",])),
        "eigenvector" = c(sum(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.parametric",]$p_value_two_side<0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.parametric",]),
                          sum(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.network",]$p_value_two_side <0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.network",])),
        "Alters" = c(sum(R[R$FemPhenotypeBias == F & R$test %in% "Alters.parametric",]$p_value_two_side<0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Alters.parametric",]),
                     sum(R[R$FemPhenotypeBias == F & R$test %in% "Alters.network",]$p_value_two_side <0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Alters.network",]))

      )

      d4 = data.frame(
        "Error Type" = rep("False positives rates", 2),
        "Biases" = rep(TRUE, 2),
        "GI" = rep(TRUE, 2),
        "two.side" = rep(TRUE, 2),
        "approches" = c("Parametric", "Nertwork permutations"),
        "strength" = c(sum(R[R$FemPhenotypeBias == F & R$test %in% "Strength.parametric.corrected",]$p_value_two_side<0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Strength.parametric.corrected",]),
                       sum(R[R$FemPhenotypeBias == F & R$test %in% "Strength.network.corrected",]$p_value_two_side <0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Strength.network.corrected",])),
        "eigenvector" = c(sum(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.parametric.corrected",]$p_value_two_side<0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.parametric.corrected",]),
                          sum(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.network.corrected",]$p_value_two_side <0.05)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Eigen.network.corrected",])),
        "Alters" = c(sum(R[R$FemPhenotypeBias == F & R$test %in% "Alters.parametric.corrected3",]$p_value_two_side<0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Alters.parametric.corrected",]),
                     sum(R[R$FemPhenotypeBias == F & R$test %in% "Alters.network.corrected",]$p_value_two_side <0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F & R$test %in% "Alters.network.corrected",]))

      )

      RESULTS = rbind(d1,d2, d3, d4)
      print(RESULTS)
    }
  }
}
all(R$coef < 0)
beepr::beep()

## Simulations without biases of observation-------------------
Mat[,2]<-rep(1,NumCombinations) # no obs bias keep constant to 1
R2 =  NULL
a=b=c=1
for (a in 1:length(FemPhenotypeBias))
{
  for(b in 1:nrow(Mat))
  {
    Result<-list()
    result = NULL
    for(c in 1:nSim)
    {
      cat("#################################################################################", '\n')
      cat("Simulation: ", b, "\n")

      df = Simulation(
        GS = Mat[b,1],ObsBias = Mat[b,2], FemSexRatio = Mat[b,3],FemPhenotypeBias = FemPhenotypeBias[a], nfocals = Mat[b,4],
        N.Perm = 100)

      df$GS = Mat[b,1]
      df$ObsBias = Mat[b,2]
      df$FemSexRatio = Mat[b,3]
      df$FemPhenotypeBias = FemPhenotypeBias[a]
      df$nfocals =Mat[b,4]
      R2 = rbind(R2, df)

      cat("#################################################################################", '\n')
      cat("Parametric true positive rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Double permutation true positive rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat('\n')

      cat("Parametric true positive rates for non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Double permutation true positive rates for  non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat('\n')

      cat("Parametric true positive rates for non GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.parametric<0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for non GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for non GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Double permutation true positive rates for  non GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat('\n')

      cat("#################################################################################", '\n')
      cat("Parametric true positive rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.parametric.corrected<0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Double permutation true positive rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == T,]$Strength.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat('\n')

      cat("Parametric true positive rates for GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.parametric.corrected<0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Double permutation true positive rates for  GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == T,]$Eigen.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat('\n')

      cat("Parametric true positive rates for GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.parametric.corrected<0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Network permutation true positive rates for GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Pre-network permutation true positive rates for GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat("Double permutation true positive rates for  GI alters: ", sum(R2[R2$FemPhenotypeBias == T,]$Alters.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]), "\n")
      cat('\n')

      cat("#################################################################################", '\n')
      cat("Parametric true negatives rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation true negatives rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation true negatives rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Double permutation true negatives rates for non GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat('\n')

      cat("Parametric true negatives rates for non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation true negatives rates for non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation true negatives rates non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Double permutation true negatives rates for  non GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat('\n')

      cat("Parametric true negatives rates for non GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.parametric<0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation true negatives rates for non GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation true negatives rates for non GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Double permutation true negatives rates for  non GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat('\n')

      cat("#################################################################################", '\n')
      cat("Parametric true negatives rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.parametric.corrected<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation true negatives rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation true negatives rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Double permutation true negatives rates for GI strength: ", sum(R2[R2$FemPhenotypeBias == F,]$Strength.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat('\n')

      cat("Parametric true negatives rates for GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.parametric.corrected<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation true negatives rates for GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation true negatives rates GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Double permutation true negatives rates for  GI eigenvector: ", sum(R2[R2$FemPhenotypeBias == F,]$Eigen.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat('\n')

      cat("Parametric true negatives rates for GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.parametric.corrected<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Network permutation true negatives rates for GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Pre-network permutation true negatives rates for GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat("Double permutation true negatives rates for  GI alters: ", sum(R2[R2$FemPhenotypeBias == F,]$Alters.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]), "\n")
      cat('\n')
    }
  }
}
beepr::beep()


#######################################
##### Results
#######################################
d1 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R[R$FemPhenotypeBias == T,]$Strength.parametric>0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                 sum(R[R$FemPhenotypeBias == T,]$Strength.network >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                 sum(R[R$FemPhenotypeBias == T,]$Strength.pre.network >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                 sum(R[R$FemPhenotypeBias == T,]$Strength.double >0.05)*100/nrow(R[R$FemPhenotypeBias == T,])),
  "eigenvector" = c(sum(R[R$FemPhenotypeBias == T,]$Eigen.parametric>0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                    sum(R[R$FemPhenotypeBias == T,]$Eigen.network >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                    sum(R[R$FemPhenotypeBias == T,]$Eigen.pre.network >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                    sum(R[R$FemPhenotypeBias == T,]$Eigen.double >0.05)*100/nrow(R[R$FemPhenotypeBias == T,])),
  "Alters" = c(sum(R[R$FemPhenotypeBias == T,]$Alters.parametric>0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T,]),
               sum(R[R$FemPhenotypeBias == T,]$Alters.network >0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T,]),
               sum(R[R$FemPhenotypeBias == T,]$Alters.pre.network >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
               sum(R[R$FemPhenotypeBias == T,]$Alters.double >0.05)*100/nrow(R[R$FemPhenotypeBias == T,])),
  "Error Type" = rep("False negatives rates", 4),

  "GI" = rep(FALSE, 4)
)

d2 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R[R$FemPhenotypeBias == T,]$Strength.parametric.corrected>0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                 sum(R[R$FemPhenotypeBias == T,]$Strength.network.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                 sum(R[R$FemPhenotypeBias == T,]$Strength.pre.network.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                 sum(R[R$FemPhenotypeBias == T,]$Strength.double.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,])),
  "eigenvector" = c(sum(R[R$FemPhenotypeBias == T,]$Eigen.parametric.corrected>0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                    sum(R[R$FemPhenotypeBias == T,]$Eigen.network.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                    sum(R[R$FemPhenotypeBias == T,]$Eigen.pre.network.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
                    sum(R[R$FemPhenotypeBias == T,]$Eigen.double.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,])),
  "Alters" = c(sum(R[R$FemPhenotypeBias == T,]$Alters.parametric.corrected>0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T,]),
               sum(R[R$FemPhenotypeBias == T,]$Alters.network.corrected >0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == T,]),
               sum(R[R$FemPhenotypeBias == T,]$Alters.pre.network.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,]),
               sum(R[R$FemPhenotypeBias == T,]$Alters.double.corrected >0.05)*100/nrow(R[R$FemPhenotypeBias == T,])),
  "Error Type" = rep("False negatives rates", 4),
  "Biases" = rep(TRUE, 4),
  "GI" = rep(TRUE, 4)
)


d3 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R[R$FemPhenotypeBias == F,]$Strength.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                 sum(R[R$FemPhenotypeBias == F,]$Strength.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                 sum(R[R$FemPhenotypeBias == F,]$Strength.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                 sum(R[R$FemPhenotypeBias == F,]$Strength.double <0.05)*100/nrow(R[R$FemPhenotypeBias == F,])),
  "eigenvector" = c(sum(R[R$FemPhenotypeBias == F,]$Eigen.parametric<0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                    sum(R[R$FemPhenotypeBias == F,]$Eigen.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                    sum(R[R$FemPhenotypeBias == F,]$Eigen.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                    sum(R[R$FemPhenotypeBias == F,]$Eigen.double <0.05)*100/nrow(R[R$FemPhenotypeBias == F,])),
  "Alters" = c(sum(R[R$FemPhenotypeBias == F,]$Alters.parametric<0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F,]),
               sum(R[R$FemPhenotypeBias == F,]$Alters.network <0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F,]),
               sum(R[R$FemPhenotypeBias == F,]$Alters.pre.network <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
               sum(R[R$FemPhenotypeBias == F,]$Alters.double <0.05)*100/nrow(R[R$FemPhenotypeBias == F,])),
  "Error Type" = rep("False positives rates", 4),
  "Biases" = rep(TRUE, 4),
  "GI" = rep(FALSE, 4)
)


d4 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R[R$FemPhenotypeBias == F,]$Strength.parametric.corrected<0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                 sum(R[R$FemPhenotypeBias == F,]$Strength.network.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                 sum(R[R$FemPhenotypeBias == F,]$Strength.pre.network.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                 sum(R[R$FemPhenotypeBias == F,]$Strength.double.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,])),
  "eigenvector" = c(sum(R[R$FemPhenotypeBias == F,]$Eigen.parametric.corrected<0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                    sum(R[R$FemPhenotypeBias == F,]$Eigen.network.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                    sum(R[R$FemPhenotypeBias == F,]$Eigen.pre.network.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
                    sum(R[R$FemPhenotypeBias == F,]$Eigen.double.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,])),
  "Alters" = c(sum(R[R$FemPhenotypeBias == F,]$Alters.parametric.corrected<0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F,]),
               sum(R[R$FemPhenotypeBias == F,]$Alters.network.corrected <0.05, na.rm = T)*100/nrow(R[R$FemPhenotypeBias == F,]),
               sum(R[R$FemPhenotypeBias == F,]$Alters.pre.network.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,]),
               sum(R[R$FemPhenotypeBias == F,]$Alters.double.corrected <0.05)*100/nrow(R[R$FemPhenotypeBias == F,])),
  "Error Type" = rep("False positives rates", 4),
  "Biases" = rep(TRUE, 4),
  "GI" = rep(TRUE, 4)
)

d5 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R2[R2$FemPhenotypeBias == T,]$Strength.parametric>0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                 sum(R2[R2$FemPhenotypeBias == T,]$Strength.network >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                 sum(R2[R2$FemPhenotypeBias == T,]$Strength.pre.network >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                 sum(R2[R2$FemPhenotypeBias == T,]$Strength.double >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,])),
  "eigenvector" = c(sum(R2[R2$FemPhenotypeBias == T,]$Eigen.parametric>0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                    sum(R2[R2$FemPhenotypeBias == T,]$Eigen.network >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                    sum(R2[R2$FemPhenotypeBias == T,]$Eigen.pre.network >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                    sum(R2[R2$FemPhenotypeBias == T,]$Eigen.double >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,])),
  "Alters" = c(sum(R2[R2$FemPhenotypeBias == T,]$Alters.parametric>0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
               sum(R2[R2$FemPhenotypeBias == T,]$Alters.network >0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
               sum(R2[R2$FemPhenotypeBias == T,]$Alters.pre.network >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
               sum(R2[R2$FemPhenotypeBias == T,]$Alters.double >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,])),
  "Error Type" = rep("False negatives rates", 4),
  "Biases" = rep(FALSE, 4),
  "GI" = rep(FALSE, 4)
)

d6 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R2[R2$FemPhenotypeBias == T,]$Strength.parametric.corrected>0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                 sum(R2[R2$FemPhenotypeBias == T,]$Strength.network.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                 sum(R2[R2$FemPhenotypeBias == T,]$Strength.pre.network.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                 sum(R2[R2$FemPhenotypeBias == T,]$Strength.double.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,])),
  "eigenvector" = c(sum(R2[R2$FemPhenotypeBias == T,]$Eigen.parametric.corrected>0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                    sum(R2[R2$FemPhenotypeBias == T,]$Eigen.network.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                    sum(R2[R2$FemPhenotypeBias == T,]$Eigen.pre.network.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
                    sum(R2[R2$FemPhenotypeBias == T,]$Eigen.double.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,])),
  "Alters" = c(sum(R2[R2$FemPhenotypeBias == T,]$Alters.parametric.corrected>0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
               sum(R2[R2$FemPhenotypeBias == T,]$Alters.network.corrected >0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
               sum(R2[R2$FemPhenotypeBias == T,]$Alters.pre.network.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,]),
               sum(R2[R2$FemPhenotypeBias == T,]$Alters.double.corrected >0.05)*100/nrow(R2[R2$FemPhenotypeBias == T,])),
  "Error Type" = rep("False negatives rates", 4),
  "Biases" = rep(FALSE, 4),
  "GI" = rep(TRUE, 4)
)


d7 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R2[R2$FemPhenotypeBias == F,]$Strength.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                 sum(R2[R2$FemPhenotypeBias == F,]$Strength.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                 sum(R2[R2$FemPhenotypeBias == F,]$Strength.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                 sum(R2[R2$FemPhenotypeBias == F,]$Strength.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,])),
  "eigenvector" = c(sum(R2[R2$FemPhenotypeBias == F,]$Eigen.parametric<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                    sum(R2[R2$FemPhenotypeBias == F,]$Eigen.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                    sum(R2[R2$FemPhenotypeBias == F,]$Eigen.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                    sum(R2[R2$FemPhenotypeBias == F,]$Eigen.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,])),
  "Alters" = c(sum(R2[R2$FemPhenotypeBias == F,]$Alters.parametric<0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
               sum(R2[R2$FemPhenotypeBias == F,]$Alters.network <0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
               sum(R2[R2$FemPhenotypeBias == F,]$Alters.pre.network <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
               sum(R2[R2$FemPhenotypeBias == F,]$Alters.double <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,])),
  "Error Type" = rep("False positives rates", 4),
  "Biases" = rep(FALSE, 4),
  "GI" = rep(FALSE, 4)
)


d8 = data.frame(
  "approches" = c("Parametric", "Nertwork permutations", "Pre-network permutation", "Double permutation"),
  "strength" = c(sum(R2[R2$FemPhenotypeBias == F,]$Strength.parametric.corrected<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                 sum(R2[R2$FemPhenotypeBias == F,]$Strength.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                 sum(R2[R2$FemPhenotypeBias == F,]$Strength.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                 sum(R2[R2$FemPhenotypeBias == F,]$Strength.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,])),
  "eigenvector" = c(sum(R2[R2$FemPhenotypeBias == F,]$Eigen.parametric.corrected<0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                    sum(R2[R2$FemPhenotypeBias == F,]$Eigen.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                    sum(R2[R2$FemPhenotypeBias == F,]$Eigen.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
                    sum(R2[R2$FemPhenotypeBias == F,]$Eigen.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,])),
  "Alters" = c(sum(R2[R2$FemPhenotypeBias == F,]$Alters.parametric.corrected<0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
               sum(R2[R2$FemPhenotypeBias == F,]$Alters.network.corrected <0.05, na.rm = T)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
               sum(R2[R2$FemPhenotypeBias == F,]$Alters.pre.network.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,]),
               sum(R2[R2$FemPhenotypeBias == F,]$Alters.double.corrected <0.05)*100/nrow(R2[R2$FemPhenotypeBias == F,])),
  "Error Type" = rep("False positives rates", 4),
  "Biases" = rep(FALSE, 4),
  "GI" = rep(TRUE, 4)
)

RESULTS = rbind(d1, d2, d3, d4, d5, d6, d7, d8)

write.csv(RESULTS, file = "results simulation1.csv")
