library(ANTs)
library(ggplot2)
library(ggpubr)
########################################################
# SRI for directed behaviours----------------------
directed.sri <- function(df, scan, actor = 'Actor', receiver = 'Receiver', weigth = "weigth",
                         method = 'sri', ynull = FALSE){
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
  Xab = df.to.mat(df, col_actor, col_receiver, col_w)
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
    mtmp = df.to.mat(tmp[[a]], col_actor, col_receiver, sym = T) 
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
      # index of interactions corrected by samplign effort
      result = (Xab/(Xab + Xba + Ya + Yb + Yab + Ynull))/(Xab + Xba + Ya  + Yb )
    }else{
      # index of interactions corrected by samplign effort
      result = (Xab/(Xab + Xba + Ya + Yb + Yab ))/(Xab + Xba + Ya  + Yb )
    }
    
    result[is.na(result)] = 0
    return(result)
  }else{
    "Currently not available."
  }
  return(result)
}

# Create preferential interactions ---------------------------------------------------------------
#' @param ids a numeric vector of individuals in the network
#' @param density network density to reach
#' @param non.alters.prob probability of creation of links to assign between non-alters
#' @param Ngroups number of subgroups within the population. The idea is to create subgroups with higher probability of occurence as the index of interaction accounts for the presence but not for interactions
#' @param groups.size.scan Mean group size per scan based on a normal distribution of standard deviance of 2
pref.attch <- function (ids  = 1:50, density = 0.2, non.alters.prob = 0.1, Ngroups= 4,
                        groups.size.scan = 6, return.grp = FALSE){
  #if(density > 1){stop("Argument denstiy need to be higer than 0 an lower than 1")}
  if(non.alters.prob > 1){stop("Argument denstiy need to be higer than 0 an lower than 1")}
  # Creating groups within the population of individuals that will be observed together but not neceserlly interacting together.
  ids.group = sample(1:Ngroups, length(ids), replace = T)
  
  m = matrix(0, length(ids), length(ids)) # Matrix size
  colnames(m) = rownames(m) = ids
  
  density2 = 0 # density to reach
  a = 1
  result = NULL
  
  while (density2 < density) {
    
    scan = a
    # Pick a group
    grp.sample = sample(1:Ngroups, 1)
    
    # Individuals in the group
    ids.group.2 = ids[ids.group %in% grp.sample]
    
    # Number of individuals observe in this scan
    N = round(rnorm(n = 1, mean = groups.size.scan, sd = 2))
    
    if(length(ids.group.2) == 1){
      egos = ids.group.2
      alters = sample(ids[!ids %in% egos], 1)
      
      if(N > 1){
        extra.egos = sample(ids[!ids %in% c(egos, alters)], 1)
      }
    }else{
      # If N is higher than total individuals in the selected group
      if(N > length(ids.group.2)){
        # divide group by 2 (givers and receivers)
        egos =  sample(ids.group.2, length(ids.group.2)/2)
        
        alters = ids.group.2[!ids.group.2 %in% egos]
        
        # Pick in the population the amount of individuals to fill the group. This will, link between the different groups in the populations and will create 'noise' in interactions
        # Those egos will receive rather than give to simulate kind of exclusion from the groups within the population
        extra.egos = sample(ids[!ids.group %in% grp.sample], ceiling((N/2) - length(egos)))
        
      }else{
        extra.egos = NULL
        # divide group by 2 (givers and receivers)
        egos =  sample(ids.group.2, length(ids.group.2)/2)
        alters = ids.group.2[!ids.group.2 %in% egos]
      }
    }
    
    
    
    # For loop for egos 
    for (b in 1:length(egos)) {
      # Select one individual
      ego = egos[b]
      
      # Check ego alters interactions 
      emmited.agression = m[ego, alters]
      
      if(length(emmited.agression) == 1){
        names(emmited.agression) = alters
      }
      
      # If no interactions exist between individuals, we sample an individual within the group
      if(all(emmited.agression == 0)){
        non.alters = as.numeric(names(emmited.agression))
        if (length(non.alters) > 1 ) { alter = sample(non.alters, 1) } else { alter = non.alters }
      }else{
        # Evaluate if ego will interact with previous interacted individual or with new one
        test = sample(0:1, 1, prob = c(non.alters.prob, (1-non.alters.prob)))
        # If ego interact with alters
        if(test == 1){
          past.alters = as.numeric(names(emmited.agression)[which(emmited.agression != 0)])
          alter = if (length(past.alters) > 1 ) { alter = sample(past.alters, 1) } else { alter = past.alters }
        }else{
          # all interaction within the group may happen, if so we pick randomly and individual within the population
          if(all(emmited.agression != 0)){
            alter = sample(ids[!ids %in% ego],1)
          }else{
            non.alters = as.numeric(names(emmited.agression)[which(emmited.agression == 0)])
            alter = if (length(non.alters) > 1 ) { alter = sample(non.alters, 1) } else { alter = non.alters }
          }
        }
      }
      if(ego == alter){stop(print("same ego later"))}
      # Adding interactions in the matrix for next interaction creations
      m[ego, alter] = m[ego, alter] + 1
      # Updating data frame of interactions
      result = rbind(result, data.frame(ego,alter, scan))
    }
    
    if(!is.null(extra.egos)){
      # For loop for extra.egos
      for (b in 1:length(extra.egos)) {
        # Select one individual
        ego = extra.egos[b]
        
        # Check ego alters interactions 
        received.agression = m[alters, ego]
        
        if(length(received.agression) == 1){
          names(received.agression) = alters
        }
        
        # If no interactions exist between individuals, we sample an individual within the group
        if(all(received.agression == 0)){
          non.alters = as.numeric(names(received.agression))
          if (length(non.alters) > 1 ) { alter = sample(non.alters, 1) } else { alter = non.alters }
        }else{
          # Evaluate if ego will interact with previous interacted individual or with new one
          test = sample(0:1, 1, prob = c(non.alters.prob, (1-non.alters.prob)))
          # If ego interact with previous interacted individual
          if(test == 1){
            past.alters = as.numeric(names(received.agression)[which(received.agression != 0)])
            if (length(past.alters) > 1 ) { alter = sample(past.alters, 1) } else { alter = past.alters }
          }else{
            if(all(received.agression != 0)){
              alter = sample(ids[!ids %in% ego],1)
            }else{
              non.alters = as.numeric(names(received.agression)[which(received.agression == 0)])
              if (length(non.alters) > 1 ) { alter = sample(non.alters, 1) } else { alter = non.alters }
            }
          }
        }
        # Adding interactions in the matrix for next interaction creations
        m[alter, ego]  = m[alter, ego] + 1
        # Updating data frame of interactions
        result = rbind(result, data.frame(ego,alter, scan))
      }
    }
    
    # Computing density
    density2 = met.density(m)
    
    a = a +1
    #cat("link creation: ", a, "density: ", density2, "\r")
  }
  if(return.grp){
    return(list("data" = result, "grp" = ids.group))
  }
  return(result)
}
##################################################

# The simulation          1) creates a dataset, 
#                         2) computes individuals' outstrength and corrects it by the time of observation, 
#                         3) creates a variable related to individuals' outstrength,
#                         4) randomizes the same variable to  create a non-relationship between individuals' variable and outstrength
#                         5) Performs node label permutation tests (you will find the commented permutation approach)
simulation <- function(ids =1:50, density = 0.4, non.alters.prob = 0.1, ynull = F, nperm = 1000, bias = 10){
  # Data creation----------
  ## Creating data set of non random associations
  dat2 = pref.attch(ids =ids, density = density, non.alters.prob = non.alters.prob)
  colnames(dat2) = c("ego", "alter", "scan")
  dat2$weigth = 1
  dat2$loc = "loc"
  ## Time of observation per dinividuals
  ids =unique(c(dat2$ego, dat2$alter))
  nobs = data.frame("id" = ids)
  nobs$nobs = NA
  for (b in 1:length(ids)) {
    nobs[nobs$id %in% ids[b],]$nobs = length(unique(dat2[dat2$ego %in% ids[b] | dat2$alter %in% ids[b],]$scan))
  }
  #perms = perm.ds.directed(df = dat2, scan= "scan", ctrlf= "loc", actor = 1, receiver = 2, nperm = nperm)
  ## Simulating relationship with explanatory variable and realizing Linear model--------------
  # Compute index of interactions
  m = directed.sri(dat2, scan= "scan",  actor = 1, receiver = 2, ynull = F)
  d2 = df.create(m)
  d2$outstrength = met.outstrength(m, dfid = "id")
  d2$degree = met.degree(m, sym = FALSE, dfid = "id")
  d2$eigen = met.eigen(m, sym = FALSE, dfid = "id")
  d2 = merge(d2, nobs, by = "id", all = T)
  d2$degree = d2$degree / d2$nobs

  # Creation of a variable that is related to individuals outstrength 
  d2 = d2[order(d2$outstrength),]
  trait = rnorm(nrow(d2),0,2)
  d2$trait = trait[order(trait)]
  d2$trait = d2$trait-min(d2$trait) + 1
  d2$trait.rand = sample(d2$trait)
  test = T
  
  # Randomize the same variable to  create a non relationship between individuals variable and outstrength
  while(test){
    # testing non random network with random relationship with y variable
    d2$trait.rand = sample(d2$trait) 
    
    # We want a non significant random association
    sig = summary(lm(outstrength ~ trait.rand, data = d2))$coefficients[2,4] > 0.05
    # We want positive relationship to facilitate p-value computation
    coef.sign = summary(lm(outstrength ~ trait.rand, data = d2))$coefficients[2,1] > 0
    
    if(all(c(coef.sign, sig) == T)){test = FALSE}
  }

  # Incorporating biases-------------
  bias = bias /1:nrow(d2)
  bias = rev(bias)
  bias
  to.remove = round((bias*d2$nobs)/100)
  to.remove
  if(all(to.remove == 0)){new.dat = dat2}else{
    new.dat = NULL
    for (a in 1:nrow(d2)) {
      if(to.remove[a] == 0){next()}
      ego = d2$id[a]
      tmp = dat2[dat2$ego %in% ego | dat2$alter %in% ego, ]
      tmp = tmp[-sample(1:nrow(tmp), to.remove[a], replace = FALSE),]
      new.dat = rbind(new.dat, tmp)
    }
  }

  # Compute index of interactions
  m2 = directed.sri(new.dat, scan= "scan",  actor = 1, receiver = 2, ynull = FALSE)
  tmp = df.create(m2)
  tmp = met.outstrength(m2, df=tmp, dfid = "id")
  tmp = met.degree(m2, df=tmp, sym = FALSE, dfid = "id")
  tmp = met.eigen(m2, df=tmp, sym = FALSE, dfid = "id")
  tmp = tmp[,-ncol(tmp)]# bouble column insert by ANTs function (ANTs)
  colnames(tmp) = paste(colnames(tmp),".bias", sep = "")
  colnames(tmp)[1] = "id"
  d2 = merge(d2, tmp, by = "id")
  d2$degree.bias = d2$degree.bias / d2$nobs
  summary(lm(outstrength.bias ~ trait, data = d2))

  p0 = ggplot(d2, aes(x = outstrength, y = trait))+geom_point()
  p1 = ggplot(d2, aes(x = outstrength.bias, y = trait))+geom_point()
  p2 = ggplot(d2, aes(x = outstrength.bias, y = outstrength))+
    geom_point(aes(x = outstrength.bias, y = outstrength))
  print(ggarrange(p0, p1,p2, nrow = 1, ncol = 3))
  # Node label ---------------

  
  r = c(coefficients(lm(outstrength.bias ~ trait, data = d2))[2], 
        coefficients(lm(outstrength.bias ~ trait.rand, data = d2))[2],
        
        coefficients(lm(degree.bias ~ trait, data = d2))[2], 
        coefficients(lm(degree.bias ~ trait.rand, data = d2))[2],
        
        coefficients(lm(eigen.bias ~ trait, data = d2))[2], 
        coefficients(lm(eigen.bias ~ trait.rand, data = d2))[2])

  for (b in 1:nperm) {
    r1 = coefficients(lm(outstrength.bias ~ sample(trait), data = d2))[2]
    r2 = coefficients(lm(outstrength.bias ~ sample(trait.rand), data = d2))[2]
    
    r3 = coefficients(lm(degree.bias ~ sample(trait), data = d2))[2]
    r4 = coefficients(lm(degree.bias ~ sample(trait.rand), data = d2))[2]
    
    r5 = coefficients(lm(eigen.bias ~ sample(trait), data = d2))[2]
    r6 = coefficients(lm(eigen.bias ~ sample(trait.rand), data = d2))[2]
    
    r = rbind(r, c(r1,r2,r3,r4,r5,r6))
  }
  r = as.data.frame(r)
  r[,7] = summary(lm(outstrength.bias ~ trait, data = d2))$coefficients[2,4]
  r[,8] = summary(lm(outstrength.bias ~ trait.rand, data = d2))$coefficients[2,4]
  r[,9] = summary(lm(degree.bias ~ trait, data = d2))$coefficients[2,4]
  r[,10] = summary(lm(degree.bias ~ trait.rand, data = d2))$coefficients[2,4]
  r[,11] = summary(lm(eigen.bias ~ trait, data = d2))$coefficients[2,4]
  r[,12] = summary(lm(eigen.bias ~ trait.rand, data = d2))$coefficients[2,4]
  rownames(r) = NULL
  colnames(r) = c("Outstrength.trait", "Outstrength.trait.rand", 
                  "Degree.trait", "Degree.trait.rand", 
                  "Eigen.trait", "Eigen.trait.rand", 
                  "Outstrengthparametric.p.non.rand", "Outstrengthparametric.p.rand",
                  "Degree.parametric.p.non.rand", "Degree.parametric.p.rand",
                  "Eigen.parametric.p.non.rand", "Eigen.parametric.p.rand")
  #R0$perm = 1:nrow(R0)
  r$perm = 1:nrow(r)
  return(r)
}

# Latin hypercube sampling--------------------------------------
library(lhs)
NumCombinations<-500

VariablesToSample<-4
VarNames<-c("GroupSize",    ## Range 10-100
            "density",  ## Range 0.2-0.80
            "non.alters.prob",## Range 0.1-0.3
            "biases")## Range 1-20

LHS<-randomLHS(NumCombinations,VariablesToSample)
Mat<-matrix(NA,nrow=NumCombinations,ncol=VariablesToSample)
Mat[,1]<-round((30 + (LHS[,1]*(100-10))),0)
Mat[,2]<-round(0.2 + (LHS[,2]*(0.8-0.2)),2)
Mat[,3]<-round(0.2 + (LHS[,3]*(0.25-0.1)),2)
Mat[,4]<-round(0.2 + (LHS[,4]*(40-1)),2)
#Mat[,4] = ifelse(Mat[,4]>1, 1, Mat[,4])
#Mat[,5]<-round(0.2 + (LHS[,4]*(40-1)),2)

a = 1
result1 = result2 = result3 =  NULL
for (a in a:nrow(Mat)) {
  cat("Simulation ", a,
      ", N individuals = ", Mat[a,1],
       ", biases = ", Mat[a,4],"%, ",
      "preferential attachment = ", (1 - Mat[a,3])*100,"%","\n")
  tmp =  simulation(ids = 1:Mat[a,1], density = Mat[a,2], non.alters.prob = Mat[a,3], ynull = F, nperm = 10000, bias = Mat[a,4])

  tmp2 = ANTs:::stat.p(tmp$Outstrength.trait)
  tmp3 = ANTs:::stat.p(tmp$Outstrength.trait.rand)
  tmp4 = as.data.frame(rbind(tmp2, tmp3))
  tmp4$type = c("non random", "random")
  tmp4$p = c(tmp$Outstrength.parametric.p.non.rand[1], tmp$Outstrength.parametric.p.rand[1])
  result1 = rbind(result1,tmp4)
  
  tmp4 = ANTs:::stat.p(tmp$Degree.trait)
  tmp5 = ANTs:::stat.p(tmp$Degree.trait.rand)
  tmp6 = as.data.frame(rbind(tmp4, tmp5))
  tmp6$type = c("non random", "random")
  tmp6$p = c(tmp$Degree.parametric.p.non.rand[1], tmp$Degree.parametric.p.rand[1])
  result2 = rbind(result2,tmp6)
  
  tmp6 = ANTs:::stat.p(tmp$Eigen.trait)
  tmp7 = ANTs:::stat.p(tmp$Eigen.trait.rand)
  tmp8 = as.data.frame(rbind(tmp6, tmp7))
  tmp8$type = c("non random", "random")
  tmp8$p = c(tmp$Eigen.parametric.p.non.rand[1], tmp$Eigen.parametric.p.rand[1])
  result3 = rbind(result3,tmp8)

  
  
  cat("Rates of false positives for outstrength with GI approach: ",
      (nrow(result1[result1$type %in% "random" & result1$`p-value_one_side` < 0.05, ])*100)/nrow(result1[result1$type %in% "random",]), "\n")
  cat("Rates of false negatives for outstrength with GI approach: ",
      (nrow(result1[result1$type %in% "non random" & result1$`p-value_one_side` > 0.05, ])*100)/nrow(result1[result1$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "random" & result2$`p-value_one_side` < 0.05, ])*100)/nrow(result2[result2$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "non random" & result2$`p-value_one_side` > 0.05, ])*100)/nrow(result2[result2$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for eigenvector with GI approach: ",
      (nrow(result3[result3$type %in% "random" & result3$`p-value_one_side` < 0.05, ])*100)/nrow(result3[result3$type %in% "random",]), "\n")
  cat("Rates of false negatives for eigenvector with GI approach: ",
      (nrow(result3[result3$type %in% "non random" & result3$`p-value_one_side` > 0.05, ])*100)/nrow(result3[result3$type %in% "non random",]), "\n")
  
  
  cat("\n")
  # Rates of false negatives
  cat("Parametric rates of false positives for outstrength with GI approach: ",
      (nrow(result1[result1$type %in% "random" & result1$p < 0.05, ])*100)/nrow(result1[result1$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for outstrength with GI approach: ",
      (nrow(result1[result1$type %in% "non random" & result1$p > 0.05, ])*100)/nrow(result1[result1$type %in% "non random",]), "\n")
  
  cat("Parametric rates of false positives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "random" & result2$p < 0.05, ])*100)/nrow(result2[result2$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "non random" & result2$p > 0.05, ])*100)/nrow(result2[result2$type %in% "non random",]), "\n")
  
  cat("Parametric rates of false positives for eigenvector with GI approach: ",
      (nrow(result3[result3$type %in% "random" & result3$p < 0.05, ])*100)/nrow(result3[result3$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for eigenvector with GI approach: ",
      (nrow(result3[result3$type %in% "non random" & result3$p > 0.05, ])*100)/nrow(result3[result3$type %in% "non random",]), "\n")
  cat("\n")
}

# No biases ---------------------
Mat[,4] = 0
a = 1
result4 = result5 = result6 =  NULL
for (a in a:nrow(Mat)) {
  cat("Simulation ", a,
      ", N individuals = ", Mat[a,1],
      ", biases = ", Mat[a,4],"%, ",
      "preferential attachment = ", (1 - Mat[a,3])*100,"%","\n")
  tmp =  simulation(ids = 1:Mat[a,1], density = Mat[a,2], non.alters.prob = Mat[a,3], ynull = F, nperm = 10000, bias = Mat[a,4])
  
  tmp2 = ANTs:::stat.p(tmp$Outstrength.trait)
  tmp3 = ANTs:::stat.p(tmp$Outstrength.trait.rand)
  tmp4 = as.data.frame(rbind(tmp2, tmp3))
  tmp4$type = c("non random", "random")
  tmp4$p = c(tmp$Outstrength.parametric.p.non.rand[1], tmp$Outstrength.parametric.p.rand[1])
  result4 = rbind(result4,tmp4)
  
  tmp4 = ANTs:::stat.p(tmp$Degree.trait)
  tmp5 = ANTs:::stat.p(tmp$Degree.trait.rand)
  tmp6 = as.data.frame(rbind(tmp4, tmp5))
  tmp6$type = c("non random", "random")
  tmp6$p = c(tmp$Degree.parametric.p.non.rand[1], tmp$Degree.parametric.p.rand[1])
  result5 = rbind(result5,tmp6)
  
  tmp6 = ANTs:::stat.p(tmp$Eigen.trait)
  tmp7 = ANTs:::stat.p(tmp$Eigen.trait.rand)
  tmp8 = as.data.frame(rbind(tmp6, tmp7))
  tmp8$type = c("non random", "random")
  tmp8$p = c(tmp$Eigen.parametric.p.non.rand[1], tmp$Eigen.parametric.p.rand[1])
  result6 = rbind(result6,tmp8)
  
  
  
  cat("Rates of false positives for outstrength with GI approach: ",
      (nrow(result4[result4$type %in% "random" & result4$`p-value_one_side` < 0.05, ])*100)/nrow(result4[result4$type %in% "random",]), "\n")
  cat("Rates of false negatives for outstrength with GI approach: ",
      (nrow(result4[result4$type %in% "non random" & result4$`p-value_one_side` > 0.05, ])*100)/nrow(result4[result4$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for outdegree with GI approach: ",
      (nrow(result5[result5$type %in% "random" & result5$`p-value_one_side` < 0.05, ])*100)/nrow(result5[result5$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree with GI approach: ",
      (nrow(result5[result5$type %in% "non random" & result5$`p-value_one_side` > 0.05, ])*100)/nrow(result5[result5$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for eigenvector with GI approach: ",
      (nrow(result6[result6$type %in% "random" & result6$`p-value_one_side` < 0.05, ])*100)/nrow(result6[result6$type %in% "random",]), "\n")
  cat("Rates of false negatives for eigenvector with GI approach: ",
      (nrow(result6[result6$type %in% "non random" & result6$`p-value_one_side` > 0.05, ])*100)/nrow(result6[result6$type %in% "non random",]), "\n")
  
  cat("\n")
  # Rates of false negatives
  cat("Parametric rates  of false positives for outstrength with GI approach: ",
      (nrow(result4[result4$type %in% "random" & result4$p < 0.05, ])*100)/nrow(result4[result4$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for outstrength with GI approach: ",
      (nrow(result4[result4$type %in% "non random" & result4$p > 0.05, ])*100)/nrow(result4[result4$type %in% "non random",]), "\n")
  
  cat("Parametric rates  of false positives for outdegree with GI approach: ",
      (nrow(result5[result5$type %in% "random" & result5$p < 0.05, ])*100)/nrow(result5[result5$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for outdegree with GI approach: ",
      (nrow(result5[result5$type %in% "non random" & result5$p > 0.05, ])*100)/nrow(result5[result5$type %in% "non random",]), "\n")
  
  cat("Parametric rates  of false positives for eigenvector with GI approach: ",
      (nrow(result6[result6$type %in% "random" & result6$p < 0.05, ])*100)/nrow(result6[result6$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for eigenvector with GI approach: ",
      (nrow(result6[result6$type %in% "non random" & result6$p > 0.05, ])*100)/nrow(result6[result6$type %in% "non random",]), "\n")
  cat("\n")
}