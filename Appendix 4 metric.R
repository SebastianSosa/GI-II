library(ANTs)
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
pref.attch <- function (ids  = 50, OBS = 1000, non.alters.prob = 0.1, Ngroups= 4,
                        groups.size.scan = 6, return.grp = FALSE, p2 = 0.2){
  #if(density > 1){stop("Argument denstiy need to be higer than 0 an lower than 1")}
  if(non.alters.prob > 1){stop("Argument denstiy need to be higer than 0 an lower than 1")}
  # Creating groups within the population of individuals that will be observed together but not neceserlly interacting together.
  ids = 1:ids
  ids.group = sample(1:Ngroups, length(ids), replace = T)

  m = matrix(0, length(ids), length(ids)) # Matrix size
  colnames(m) = rownames(m) = as.character(ids)

  density2 = 0 # density to reach
  a = 1
  result = NULL

  while (a < OBS) {
    scan = a
    # Pick a group
    grp.sample = sample(1:Ngroups, 1)

    # Individuals in the group
    ids.group.2 = ids[ids.group %in% grp.sample]

    if(length(ids.group.2) == 1){
      # test if individuals from outside the group will be in the scan
      extra.egos = FALSE
      extra.egos = sample(c(FALSE, TRUE), 1, prob = c((1-p2),p2))
      if(all(extra.egos != FALSE) == TRUE){
        #cat("Grp selected of length 1", "\n")
        ego = ids.group.2
        alter = sample(ids[!ids.group %in% grp.sample], 1)
        if(ego == alter){stop(print("same ego later"))}
        m[alter, ego]  = m[alter, ego] + 1
        # Updating data frame of interactions
        result = rbind(result, data.frame(ego,alter, scan))
        next()
      }else{
        # If no extra-ego then no associations observed
        next()
      }
    }

    # Number of individuals observe in this scan
    N = rpois(1,groups.size.scan)
    while(N < 2 ){
      N = rpois(1,groups.size.scan)
    }

    # test if individuals from outside the group will be in the scan
    extra.egos = FALSE
    extra.egos = sample(c(FALSE, TRUE), 1, prob = c((1-p2),p2))
    if(extra.egos){
      # If picked group size is lower than N then extra ego will be added at the end
      if(N > length(ids.group.2)){
        # divide group by 2 (givers and receivers)
        egos =  sample(ids.group.2, length(ids.group.2)/2)

        alters = ids.group.2[!ids.group.2 %in% egos]

        # Pick in the population the amount of individuals to fill the group. This will, link between the different groups in the populations and will create 'noise' in interactions
        # Those egos will receive rather than give to simulate kind of exclusion from the groups within the population
        extra.egos = sample(ids[!ids.group %in% grp.sample], ceiling((N/2) - length(egos)))

      }else{
        # If picked group size is higher than N then an extra ego is picked at the begingin then egos from same group are picked
        extra.egos = sample(ids[!ids.group %in% grp.sample],1)
        if(N > 2){
          # divide group by 2 (givers and receivers)
          egos =  sample(ids.group.2, (N-1)/2)
          alters = ids.group.2[!ids.group.2 %in% egos]
        }else{
          # Else, observation is of one extra ego and one individual from the group
          #cat("Grp selected of length larger than 1", "\n")
          ego = sample(ids.group.2, 1)
          alter = sample(ids[!ids.group %in% grp.sample], 1)
          if(ego == alter){stop(print("same ego later"))}
          m[alter, ego]  = m[alter, ego] + 1
          # Updating data frame of interactions
          result = rbind(result, data.frame(ego,alter, scan))
          ego = alter = NA
          next()
        }
      }
    }else{
      # If picked group size is lower than N then N will not be achieved
      if(N > length(ids.group.2)){
        egos =  sample(ids.group.2, length(ids.group.2)/2)
        alters = ids.group.2[!ids.group.2 %in% egos]
      }else{
        egos =  sample(ids.group.2, N/2)
        alters = ids.group.2[!ids.group.2 %in% egos]
      }
    }

    # For loop for egos
    for (b in 1:length(egos)) {


      # Select one individual
      ego = egos[b]

      # Check ego alters interactions
      emmited.agression = m[ego, alters]

      if(!all(names(emmited.agression) %in% as.character(ids.group.2))){
        stop("alters not in group")
      }

      if(length(emmited.agression) == 1){
        names(emmited.agression) = alters
      }

      # If no interactions exist between individuals, we sample an individual within the group
      if(all(emmited.agression == 0)){
        #cat("no past interactions", "\n")
        non.alters = as.numeric(names(emmited.agression))
        if (length(non.alters) > 1 ) { alter = sample(non.alters, 1) } else { alter = non.alters }
      }else{
        ##cat("Past interactions present", "\n")
        # Evaluate if ego will interact with previous interacted individual or with new one
        test = sample(0:1, 1, prob = c(non.alters.prob, (1-non.alters.prob)))
        # If ego interact with alters
        if(test == 1){
          #cat("Ego interacting with alters", "\n")
          past.alters = as.numeric(names(emmited.agression)[which(emmited.agression != 0)])
          alter = if (length(past.alters) > 1 ) { alter = sample(past.alters, 1) } else { alter = past.alters }
        }else{
          #cat("Ego not interacting with alters", "\n")
          # all interaction within the group may happen, if so we pick randomly and individual within the population
          if(all(emmited.agression != 0)){
            #cat("All interactions within the group exist", "\n")
            next()
          }else{
            #cat("All interactions within the group exist", "\n")
            non.alters = as.numeric(names(emmited.agression)[which(emmited.agression == 0)])
            alter = if (length(non.alters) > 1 ) { alter = sample(non.alters, 1) } else { alter = non.alters }
          }
        }
      }
      #cat("scan: ", a, " ids in grp: ", ids.group.2, " egos: ", egos, " alters: ", alters, ' ego : ', ego, " alter: ", alter, "\n")
      if(alter == 0){stop()}
      if(ego == alter){stop(print("same ego later"))}
      ego.grp = ids.group[ego]
      alter.grp = ids.group[alter]
      if(ego.grp != alter.grp){stop("ego and alter are not from the same group")}
      # Adding interactions in the matrix for next interaction creations
      m[ego, alter] = m[ego, alter] + 1
      # Updating data frame of interactions
      if(alter == 0){stop()}
      result = rbind(result, data.frame(ego,alter, scan))
    }

    if(all(extra.egos != FALSE) == TRUE){
      #cat("Adding extra egos", "\n")
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
        #cat("scan: ", a, " ids in grp: ", ids.group.2, " egos: ", egos, " alters: ", alters, ' ego : ', ego, " alter: ", alter, "\n")
        if(ego == alter){stop(print("same ego later"))}
        if(alter == 0){stop()}
        # Adding interactions in the matrix for next interaction creations
        m[alter, ego]  = m[alter, ego] + 1
        # Updating data frame of interactions
        result = rbind(result, data.frame(ego,alter, scan))
      }
    }
    a = a +1
    ###cat("link creation: ", a, "density: ", density2, "\r")
  }
  if(return.grp){
    return(list("data" = result, "grp" = ids.group))
  }
  return(result)
}
###################################################
#' @title  Simulation 4 for outdegree
#' @description  The simulation 1) creates a dataset, 2) computes individuals' outdegree and corrects it by the time of observation, 3) creates a variable related to individuals' outstrength, 4) randomizes the same variable to  create a non-relationship between individuals' variable and outstrength, 5) Performs node label permutation tests. In addition, biases can be simulated
#' @param ids A numeric vector of individuals in the network
#' @param OBS Number of observations to reach
#' @param non.alters.prob Probability of creation of links to assign between non-alters
#' @param Ngroups Number of subgroups within the population. The idea is to create subgroups with higher probability of occurence as the index of interaction accounts for the presence but not for interactions
#' @param groups.size.scan Mean group size per scan based on a Poisson distribution.
#' @param ynull Term of directed SRI
#' @param nperm Number of permutations to perform
#' @param bias Amount of Biases to include
#' @param p2 Probability of non-subgroup member
#' @param metric Social network metric to compute
#' @binary.correction Does binary correction need to be applied
simulation <- function(ids =50, OBS = 100, non.alters.prob = 0.1, Ngroups = 4, groups.size.scan = 6, ynull = F, nperm = 1000, bias = 10,  p2 = 0.2, metric = "met.degree", binary.correction = TRUE, ...){
  # Data creation----------
  ## Creating data set of non random associations
  dat2 = pref.attch(ids =ids, OBS = OBS, non.alters.prob = non.alters.prob,Ngroups = Ngroups, groups.size.scan = groups.size.scan, p2 = p2, return.grp = FALSE)
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
  ## Simulating relationship with explanatory variable and realizing Linear model--------------
  # Compute index of interactions
  m = directed.sri(dat2, scan= "scan",  actor = 1, receiver = 2, ynull = F)
  d2 = df.create(m)
  d2$metric = do.call(metric, list(M = m))
  d2 = merge(d2, nobs, by = "id", all = T)
  if(binary.correction){
    d2$metric.corrected <- d2$metric/d2$nobs
  }else{d2$metric.corrected <- d2$metric}

  # Creation of a variable that is related to individuals outstrength ----------
  d2 = d2[order(d2$metric.corrected),]
  trait = rnorm(nrow(d2),0,2)
  d2$trait = trait[order(trait)]
  d2$trait = d2$trait-min(d2$trait) + 1
  d2$trait.rand = sample(d2$trait)
  test = T

  # Randomize the same variable to  create a non relationship between individuals variable and outstrength----------
  while(test){
    # testing non random network with random relationship with y variable
    d2$trait.rand = sample(d2$trait)

    # We want a non significant random association
    sig = summary(lm(metric.corrected ~ trait.rand, data = d2))$coefficients[2,4] > 0.05
    # We want positive relationship to facilitate p-value computation
    coef.sign = summary(lm(metric.corrected ~ trait.rand, data = d2))$coefficients[2,1] > 0

    if(all(c(coef.sign, sig) == T)){test = FALSE}
  }

  # Incorporating biases-------------
  bias = bias /1:nrow(d2)
  bias = rev(bias)
  bias
  to.remove = round((bias*d2$nobs)/100)
  to.remove
  if(all(to.remove == 0)){new.dat = dat2}else{
    new.dat = dat2
    for (a in 1:nrow(d2)) {
      if(to.remove[a] == 0){next()}
      ego = d2$id[a]
      tmp = which(dat2$ego %in% ego | dat2$alter %in% ego)
      new.dat = new.dat[-sample(tmp, to.remove[a], replace = FALSE),]
    }
  }

  # Compute index of interactions
  m2 = directed.sri(new.dat, scan= "scan",  actor = 1, receiver = 2, ynull = FALSE)
  tmp = df.create(m2)
  tmp$metric.bias = do.call(metric, list(M = m,...))

  ## Time of observation per dinividuals
  ids =unique(c(new.dat$ego, new.dat$alter))
  nobs = data.frame("id" = ids)
  nobs$nobs.bias = NA
  for (b in 1:length(ids)) {
    nobs[nobs$id %in% ids[b],]$nobs.bias = length(unique(new.dat[new.dat$ego %in% ids[b] | new.dat$alter %in% ids[b],]$scan))
  }
  tmp = merge(tmp, nobs, by = "id", all = T)
  if(binary.correction){
    tmp$metric.bias.corrected <- tmp$metric.bias/d2$nobs
  }else{tmp$metric.bias.corrected <- tmp$metric.bias}
  d2 = merge(d2, tmp, by = "id")


  # Node label ---------------

  r = c(coefficients(lm(metric.bias.corrected ~ trait, data = d2))[2],
        coefficients(lm(metric.bias.corrected ~ trait.rand, data = d2))[2])

  for (b in 1:nperm) {
    r1 = coefficients(lm(metric.bias.corrected ~ sample(trait), data = d2))[2]
    r2 = coefficients(lm(metric.bias.corrected ~ sample(trait.rand), data = d2))[2]
    r = rbind(r, c(r1,r2))
  }
  r = as.data.frame(r)
  r[,3] = summary(lm(metric.bias.corrected ~ trait, data = d2))$coefficients[2,4]
  cat(summary(lm(metric.bias.corrected ~ trait, data = d2))$coefficients[2,4], "\n")
  r[,4] = summary(lm(metric.bias.corrected ~ trait.rand, data = d2))$coefficients[2,4]
  rownames(r) = NULL
  colnames(r) = c(paste(metric,"trait", sep = "."),paste(metric,"trait.rand", sep = "."),
                  paste(metric,"parametric","trait", sep = "."),paste(metric,"parametric","trait.rand", sep = "."))


  #R0$perm = 1:nrow(R0)
  r$perm = 1:nrow(r)
  return(r)
}
##################################
# Latin hypercube sampling
###################
library(lhs)
NumCombinations<-500

VariablesToSample<-5
VarNames<-c("GroupSize",    ## Range 10-100
            "density",  ## Range 0.2-0.80
            "non.alters.prob",## Range 0.1-0.3
            "biases",
            "p2")## Range 1-20

LHS<-randomLHS(NumCombinations,VariablesToSample)
Mat<-matrix(NA,nrow=NumCombinations,ncol=VariablesToSample)
Mat[,1]<-round((30 + (LHS[,1]*(100-10))),0)
Mat[,2]<-round((100 + (LHS[,4]*(10000-1000))),0)
Mat[,3]<-round(0.2 + (LHS[,3]*(0.25-0.1)),2)
Mat[,4]<-round(0.2 + (LHS[,4]*(40-1)),2)
Mat[,5]<-round(0.2 + (LHS[,4]*(0.70-0.4)),2)


## Simulations with biases of observation-------------------
a = 1
result1 = NULL
for (a in a:nrow(Mat)) {
  cat("Simulation ", a,
      ", N individuals = ", Mat[a,1],
      ", biases = ", Mat[a,4],"%, ",
      "preferential attachment = ", (1 - Mat[a,3])*100,"%","\n")
  tmp =  simulation(ids = Mat[a,1], OBS = Mat[a,2], non.alters.prob = Mat[a,3], ynull = F, nperm = 10000, bias = Mat[a,4], p2 = Mat[a,5],
                    metric = "met.outstrength", binary.correction = FALSE)

  tmp4 = ANTs:::stat.p(tmp[,1])
  tmp5 = ANTs:::stat.p(tmp[,2])
  tmp6 = as.data.frame(rbind(tmp4, tmp5))
  tmp6$type = c("non random", "random")
  tmp6$p = c(tmp[,3][1], tmp[,4][1])
  result1 = rbind(result1,tmp6)


  cat("Rates of false positives for outdegree with GI approach: ",
      (nrow(result1[result1$type %in% "random" & result1$`p-value_one_side` < 0.05, ])*100)/nrow(result1[result1$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree with GI approach: ",
      (nrow(result1[result1$type %in% "non random" & result1$`p-value_one_side` > 0.05, ])*100)/nrow(result1[result1$type %in% "non random",]), "\n")


  cat("\n")
  # Rates of false negatives
  cat("Parametric rates of false positives for outdegree with GI approach: ",
      (nrow(result1[result1$type %in% "random" & result1$p < 0.05, ])*100)/nrow(result1[result1$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for outdegree with GI approach: ",
      (nrow(result1[result1$type %in% "non random" & result1$p > 0.05, ])*100)/nrow(result1[result1$type %in% "non random",]), "\n")

  cat("\n")
}

## Simulations without biases of observation-------------------
Mat[,4] = 0
a = 1
result2 =   NULL
for (a in a:nrow(Mat)) {
  cat("Simulation ", a,
      ", N individuals = ", Mat[a,1],
      ", biases = ", Mat[a,4],"%, ",
      "preferential attachment = ", (1 - Mat[a,3])*100,"%","\n")
  tmp =  simulation(ids = Mat[a,1], OBS = Mat[a,2], non.alters.prob = Mat[a,3], ynull = F, nperm = 10000, bias = Mat[a,4], p2 = Mat[a,5],
                    metric = "met.outstrength", binary.correction = FALSE)

  tmp4 = ANTs:::stat.p(tmp[,1])
  tmp5 = ANTs:::stat.p(tmp[,2])
  tmp6 = as.data.frame(rbind(tmp4, tmp5))
  tmp6$type = c("non random", "random")
  tmp6$p = c(tmp[,3][1], tmp[,4][1])
  result2 = rbind(result2,tmp6)

  cat("Rates of false positives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "random" & result2$`p-value_one_side` < 0.05, ])*100)/nrow(result2[result2$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "non random" & result2$`p-value_one_side` > 0.05, ])*100)/nrow(result2[result2$type %in% "non random",]), "\n")

  cat("\n")
  # Rates of false negatives


  cat("Parametric rates of false positives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "random" & result2$p < 0.05, ])*100)/nrow(result2[result2$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for outdegree with GI approach: ",
      (nrow(result2[result2$type %in% "non random" & result2$p > 0.05, ])*100)/nrow(result2[result2$type %in% "non random",]), "\n")

  cat("\n")
}
