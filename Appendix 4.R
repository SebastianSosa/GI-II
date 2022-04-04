library(ANTs)
############################################################################################################
######### Modification  1 (conceptual): SRI for directed behaviours
directed.sri <- function(df, scan, actor = 'Actor', receiver = 'Receiver', weigth = "weigth",
                         method = 'sri', ynull = FALSE, GI = T, return.GI = F){
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
      if(GI){result = (Xab/(Xab + Xba + Ya + Yb + Yab + Ynull))/(Xab + Xba + Ya  + Yb )}else{result = (Xab/(Xab + Xba + Ya + Yb + Yab + Ynull))}
      # index of interactions corrected by samplign effort
      
    }else{
      # index of interactions corrected by samplign effort
      if(GI){result = (Xab/(Xab + Xba + Ya + Yb + Yab ))/(Xab + Xba + Ya  + Yb )}else{result = (Xab/(Xab + Xba + Ya + Yb + Yab ))}
    }

    result[is.na(result)] = 0

  }else{
    "Currently not available."
  }
  
  if(return.GI){return(list("interaction" = result, "sampling effort" = (Xab + Xba + Ya  + Yb )))}else{return(result)}
  
}
############################################################################################################
#' @title Create preferential interactions
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
#' @param ids a numeric vector of individuals in the network
#' @param density network density to reach
#' @param non.alters.prob probability of creation of links to assign between non-alters
#' @param Ngroups number of subgroups within the population. The idea is to create subgroups with higher probability of occurence as the index of interaction accounts for the presence but not for interactions
#' @param groups.size.scan Mean group size per scan based on a normal distribution of standard deviance of 2
simulation <- function(ids =50, OBS = 1000, non.alters.prob = 0.1, ynull = F, nperm = 1000, bias = 10){
  # Data creation----------
  ## Creating data set of non random associations
  dat2 = pref.attch(ids =ids, OBS = OBS, non.alters.prob = non.alters.prob)
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
  m = directed.sri(dat2, scan= "scan",  actor = 1, receiver = 2, ynull = F, GI = F, return.GI = T)# Direct sri without correction
  d2 = df.create(m[[1]])
  d2 = merge(d2, nobs, by = "id", all = T)
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
  
  
  d2$outstrength = met.outstrength(m[[1]], dfid = "id")
  d2$degree = met.outdegree(m[[1]],  dfid = "id")
  d2$eigen = met.eigen(m[[1]], sym = F, out = T)
  
  corrected.mat = m[[1]]/m[[2]]
  corrected.mat[is.nan(corrected.mat)] = 0
  corrected.mat[is.infinite(corrected.mat)]= 0
  
  d2$outstrength.corrected = met.outstrength(corrected.mat)
  
  d2$eigen.corrected  = met.eigen(corrected.mat, sym = F, out = T)
 

  d2$degree.corrected  = d2$degree / d2$nobs

  # Creation of a variable that is related to individuals social measures--------------------
  d2 = d2[order(d2$outstrength),]
  trait = rnorm(nrow(d2),0,2)
  d2$trait.outstrength = trait[order(trait)]
  d2$trait.outstrength  = d2$trait.outstrength -min(d2$trait.outstrength ) + 1
  d2$trait.rand.outstrength  = sample(d2$trait.outstrength)
  
  d2 = d2[order(d2$degree),]
  trait = rnorm(nrow(d2),0,2)
  d2$trait.degree = trait[order(trait)]
  d2$trait.degree  = d2$trait.degree -min(d2$trait.degree ) + 1
  d2$trait.rand.degree = sample(d2$trait.degree)
  
  d2 = d2[order(d2$eigen),]
  trait = rnorm(nrow(d2),0,2)
  d2$trait.eigen = trait[order(trait)]
  d2$trait.eigen  = d2$trait.eigen -min(d2$trait.eigen) + 1
  d2$trait.rand.eigen = sample(d2$trait.eigen)


  ###### Data permutations ###############
  ############################################################################################################
  ######### Node label permutations---------------------

  r = c(coefficients(lm(outstrength ~ trait.outstrength, data = d2))[2],
        coefficients(lm(outstrength ~ trait.rand.outstrength, data = d2))[2],
        coefficients(lm(degree ~ trait.degree, data = d2))[2],
        coefficients(lm(degree ~ trait.rand.degree, data = d2))[2],
        coefficients(lm(eigen ~ trait.eigen, data = d2))[2],
        coefficients(lm(eigen ~ trait.rand.eigen, data = d2))[2],
  
        coefficients(lm(outstrength.corrected  ~ trait.outstrength, data = d2))[2],
        coefficients(lm(outstrength.corrected  ~ trait.rand.outstrength, data = d2))[2],
        coefficients(lm(degree.corrected  ~ trait.degree, data = d2))[2],
        coefficients(lm(degree.corrected ~ trait.rand.degree, data = d2))[2],
        coefficients(lm(eigen.corrected ~ trait.eigen, data = d2))[2],
        coefficients(lm(eigen.corrected ~ trait.rand.eigen, data = d2))[2]
        )
  

  for (b in 1:nperm) {
    r1 = coefficients(lm(outstrength ~ sample(trait.outstrength), data = d2))[2]
    r2 = coefficients(lm(outstrength ~ sample(trait.rand.outstrength), data = d2))[2]
    r3 = coefficients(lm(degree ~ sample(trait.degree), data = d2))[2]
    r4 = coefficients(lm(degree ~ sample(trait.rand.degree), data = d2))[2]
    r5 = coefficients(lm(eigen ~ sample(trait.eigen), data = d2))[2]
    r6 = coefficients(lm(eigen ~ sample(trait.rand.eigen), data = d2))[2]

    r7 = coefficients(lm(outstrength.corrected ~ sample(trait.outstrength), data = d2))[2]
    r8 = coefficients(lm(outstrength.corrected ~ sample(trait.rand.outstrength), data = d2))[2]
    r9 = coefficients(lm(degree.corrected ~ sample(trait.degree), data = d2))[2]
    r10 = coefficients(lm(degree.corrected ~ sample(trait.rand.degree), data = d2))[2]
    r11 = coefficients(lm(eigen.corrected ~ sample(trait.eigen), data = d2))[2]
    r12 = coefficients(lm(eigen.corrected ~ sample(trait.rand.eigen), data = d2))[2]    

    r = rbind(r, c(r1,r2,r3,r4,r5,r6, r7, r8, r9, r10, r11, r12))
  }
  r = as.data.frame(r)
  
  ############################################################################################################
  ######### Modification  1 : One-tailed parametric test
  s.degree = summary(lm(outstrength ~ trait.outstrength,data=d2))
  s.eigen = summary(lm(eigen ~ trait.eigen, data = d2))
  s.alters = summary(lm(degree ~ trait.degree, data = d2))
  
  s.degree.rand = summary(lm(outstrength ~ trait.rand.outstrength, data = d2))
  s.eigen.rand = summary(lm(eigen ~ trait.rand.eigen, data = d2))
  s.alters = summary(lm(degree ~ trait.rand.degree, data = d2))
  
  p.degree = pt(coef(s.degree)[,3], s.degree$df[2], lower = F)[2]
  p.eigen = pt(coef(s.eigen)[,3], s.eigen$df[2], lower = F)[2]
  p.alters = pt(coef(s.alters)[,3], s.alters$df[2], lower = F)[2]
  
  p.degree.rand = pt(coef(s.degree)[,3], s.degree$df[2], lower = F)[2]
  p.eigen.rand = pt(coef(s.eigen)[,3], s.eigen$df[2], lower = F)[2]
  p.alters.rand = pt(coef(s.alters)[,3], s.alters$df[2], lower = F)[2]
  
  r[,13] = p.degree
  r[,14] = p.degree.rand
  r[,15] =p.alters
  r[,16] = p.alters.rand
  r[,17] = p.eigen
  r[,18] = p.eigen.rand
  
  s.degree = summary(lm(outstrength.corrected ~ trait.outstrength,data=d2))
  s.eigen = summary(lm(eigen.corrected ~ trait.eigen, data = d2))
  s.alters = summary(lm(degree.corrected ~ trait.degree, data = d2))
  
  s.degree.rand = summary(lm(outstrength ~ trait.rand.outstrength, data = d2))
  s.eigen.rand = summary(lm(eigen ~ trait.rand.eigen, data = d2))
  s.alters = summary(lm(degree ~ trait.rand.degree, data = d2))
  
  p.degree = pt(coef(s.degree)[,3], s.degree$df[2], lower = F)[2]
  p.eigen = pt(coef(s.eigen)[,3], s.eigen$df[2], lower = F)[2]
  p.alters = pt(coef(s.alters)[,3], s.alters$df[2], lower = F)[2]
  
  p.degree.rand = pt(coef(s.degree)[,3], s.degree$df[2], lower = F)[2]
  p.eigen.rand = pt(coef(s.eigen)[,3], s.eigen$df[2], lower = F)[2]
  p.alters.rand = pt(coef(s.alters)[,3], s.alters$df[2], lower = F)[2]
  
  r[,19] = p.degree
  r[,20] = p.degree.rand
  r[,21] =p.alters
  r[,22] = p.alters.rand
  r[,23] = p.eigen
  r[,24] = p.eigen.rand
  
  rownames(r) = NULL
  colnames(r) = c("Outstrength.trait", "Outstrength.trait.rand",
                  "Degree.trait", "Degree.trait.rand",
                  "Eigen.trait", "Eigen.trait.rand",
                  "Outstrength.trait.corrected", "Outstrength.trait.rand.corrected",
                  "Degree.trait.corrected", "Degree.trait.rand.corrected",
                  "Eigen.trait.corrected", "Eigen.trait.rand.corrected",
                  
                  "Outstrength.parametric.p.non.rand", "Outstrength.parametric.p.rand",
                  "Degree.parametric.p.non.rand", "Degree.parametric.p.rand",
                  "Eigen.parametric.p.non.rand", "Eigen.parametric.p.rand",
                  "Outstrength.parametric.p.non.rand.corrected", "Outstrength.parametric.p.rand.corrected",
                  "Degree.parametric.p.non.rand.corrected", "Degree.parametric.p.rand.corrected",
                  "Eigen.parametric.p.non.rand.corrected", "Eigen.parametric.p.rand.corrected")
  #R0$perm = 1:nrow(R0)
  r$perm = 1:nrow(r)
  return(r)
}

# Latin hypercube sampling--------------------------------------
library(lhs)
NumCombinations<-500
## Simulations with biases of observation-------------------
VariablesToSample<-4
VarNames<-c("GroupSize",    ## Range 30-100
            "OBS",  ## Range 100-1000
            "non.alters.prob",## Range 0.1-0.3
            "biases")## Range 1-20

LHS<-randomLHS(NumCombinations,VariablesToSample)
Mat<-matrix(NA,nrow=NumCombinations,ncol=VariablesToSample)
Mat[,1]<-round((30 + (LHS[,1]*(100-10))),0)
Mat[,2]<-round(100 + (LHS[,2]*(1000-100)),0)
Mat[,3]<-round(0.2 + (LHS[,3]*(0.25-0.1)),2)
Mat[,4]<-round(0.2 + (LHS[,4]*(40-1)),2)

a = 1
result1 = result2 = result3 =  result4 = result5 = result6 = NULL
for (a in a:nrow(Mat)) {
  cat("#################################################################################", '\n')
  cat("Simulation ", a,
      ", N individuals = ", Mat[a,1],
       ", biases = ", Mat[a,4],"%, ",
      "preferential attachment = ", (1 - Mat[a,3])*100,"%","\n")
  tmp =  simulation(ids = Mat[a,1], OBS = Mat[a,2], non.alters.prob = Mat[a,3], ynull = F, nperm = 1000, bias = Mat[a,4])

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
  
  tmp2 = ANTs:::stat.p(tmp$Outstrength.trait.corrected)
  tmp3 = ANTs:::stat.p(tmp$Outstrength.trait.rand.corrected)
  tmp4 = as.data.frame(rbind(tmp2, tmp3))
  tmp4$type = c("non random", "random")
  tmp4$p = c(tmp$Outstrength.parametric.p.non.rand.corrected[1], tmp$Outstrength.parametric.p.rand.corrected[1])
  result4 = rbind(result4,tmp4)
  
  tmp4 = ANTs:::stat.p(tmp$Degree.trait.corrected)
  tmp5 = ANTs:::stat.p(tmp$Degree.trait.rand.corrected)
  tmp6 = as.data.frame(rbind(tmp4, tmp5))
  tmp6$type = c("non random", "random")
  tmp6$p = c(tmp$Degree.parametric.p.non.rand.corrected[1], tmp$Degree.parametric.p.rand.corrected[1])
  result5 = rbind(result5,tmp6)
  
  tmp6 = ANTs:::stat.p(tmp$Eigen.trait.corrected)
  tmp7 = ANTs:::stat.p(tmp$Eigen.trait.rand.corrected)
  tmp8 = as.data.frame(rbind(tmp6, tmp7))
  tmp8$type = c("non random", "random")
  tmp8$p = c(tmp$Eigen.parametric.p.non.rand.corrected[1], tmp$Eigen.parametric.p.rand.corrected[1])
  result6 = rbind(result6,tmp8)
  
  cat("#################################################################################", '\n')
  cat("Without GI", '\n')
  cat("Rates of false positives for outstrength approach: ",
      (nrow(result1[result1$type %in% "random" & result1$`p-value_left_side` < 0.05, ])*100)/nrow(result1[result1$type %in% "random",]), "\n")
  cat("Rates of false negatives for outstrength approach: ",
      (nrow(result1[result1$type %in% "non random" & result1$`p-value_left_side` > 0.05, ])*100)/nrow(result1[result1$type %in% "non random",]), "\n")

  cat("Rates of false positives for outdegree approach: ",
      (nrow(result2[result2$type %in% "random" & result2$`p-value_left_side` < 0.05, ])*100)/nrow(result2[result2$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree approach: ",
      (nrow(result2[result2$type %in% "non random" & result2$`p-value_left_side` > 0.05, ])*100)/nrow(result2[result2$type %in% "non random",]), "\n")

  cat("Rates of false positives for eigenvector approach: ",
      (nrow(result3[result3$type %in% "random" & result3$`p-value_left_side` < 0.05, ])*100)/nrow(result3[result3$type %in% "random",]), "\n")
  cat("Rates of false negatives for eigenvector approach: ",
      (nrow(result3[result3$type %in% "non random" & result3$`p-value_left_side` > 0.05, ])*100)/nrow(result3[result3$type %in% "non random",]), "\n")

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
  
  cat("#################################################################################", '\n')
  cat("With GI", '\n')
  cat("Rates of false positives for outstrength: ",
      (nrow(result4[result4$type %in% "random" & result4$`p-value_left_side` < 0.05, ])*100)/nrow(result4[result4$type %in% "random",]), "\n")
  cat("Rates of false negatives for outstrength approach: ",
      (nrow(result4[result4$type %in% "non random" & result4$`p-value_left_side` > 0.05, ])*100)/nrow(result1[result4$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for outdegree: ",
      (nrow(result5[result5$type %in% "random" & result5$`p-value_left_side` < 0.05, ])*100)/nrow(result5[result5$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree : ",
      (nrow(result5[result5$type %in% "non random" & result5$`p-value_left_side` > 0.05, ])*100)/nrow(result5[result5$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for eigenvector : ",
      (nrow(result6[result6$type %in% "random" & result6$`p-value_left_side` < 0.05, ])*100)/nrow(result6[result6$type %in% "random",]), "\n")
  cat("Rates of false negatives for eigenvector : ",
      (nrow(result6[result6$type %in% "non random" & result6$`p-value_left_side` > 0.05, ])*100)/nrow(result6[result6$type %in% "non random",]), "\n")
  
  # Rates of false negatives
  cat("Parametric rates of false positives for outstrength : ",
      (nrow(result4[result4$type %in% "random" & result4$p < 0.05, ])*100)/nrow(result4[result4$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for outstrength with GI approach: ",
      (nrow(result4[result4$type %in% "non random" & result4$p > 0.05, ])*100)/nrow(result4[result4$type %in% "non random",]), "\n")
  
  cat("Parametric rates of false positives for outdegree with GI : ",
      (nrow(result5[result5$type %in% "random" & result5$p < 0.05, ])*100)/nrow(result5[result5$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for outdegree with GI : ",
      (nrow(result5[result5$type %in% "non random" & result5$p > 0.05, ])*100)/nrow(result5[result5$type %in% "non random",]), "\n")
  
  cat("Parametric rates of false positives for eigenvector with GI : ",
      (nrow(result6[result6$type %in% "random" & result6$p < 0.05, ])*100)/nrow(result6[result6$type %in% "random",]), "\n")
  cat("Parametric rates of false negatives for eigenvector with GI : ",
      (nrow(result6[result6$type %in% "non random" & result6$p > 0.05, ])*100)/nrow(result6[result6$type %in% "non random",]), "\n")
  cat("\n")

}

## Simulations without biases of observation-------------------
Mat[,4] = 0
a = 1
result7 = result8 = result9 =  result10 = result11 = result12 =NULL
for (a in a:nrow(Mat)) {
  cat("Simulation ", a,
      ", N individuals = ", Mat[a,1],
      ", biases = ", Mat[a,4],"%, ",
      "preferential attachment = ", (1 - Mat[a,3])*100,"%","\n")
  tmp =  simulation(ids = Mat[a,1], OBS = Mat[a,2], non.alters.prob = Mat[a,3], ynull = F, nperm = 1000, bias = Mat[a,4])
  
  tmp2 = ANTs:::stat.p(tmp$Outstrength.trait)
  tmp3 = ANTs:::stat.p(tmp$Outstrength.trait.rand)
  tmp4 = as.data.frame(rbind(tmp2, tmp3))
  tmp4$type = c("non random", "random")
  tmp4$p = c(tmp$Outstrength.parametric.p.non.rand[1], tmp$Outstrength.parametric.p.rand[1])
  result7 = rbind(result7,tmp4)
  
  tmp4 = ANTs:::stat.p(tmp$Degree.trait)
  tmp5 = ANTs:::stat.p(tmp$Degree.trait.rand)
  tmp6 = as.data.frame(rbind(tmp4, tmp5))
  tmp6$type = c("non random", "random")
  tmp6$p = c(tmp$Degree.parametric.p.non.rand[1], tmp$Degree.parametric.p.rand[1])
  result8 = rbind(result8,tmp6)
  
  tmp6 = ANTs:::stat.p(tmp$Eigen.trait)
  tmp7 = ANTs:::stat.p(tmp$Eigen.trait.rand)
  tmp8 = as.data.frame(rbind(tmp6, tmp7))
  tmp8$type = c("non random", "random")
  tmp8$p = c(tmp$Eigen.parametric.p.non.rand[1], tmp$Eigen.parametric.p.rand[1])
  result9 = rbind(result9,tmp8)
  
  tmp2 = ANTs:::stat.p(tmp$Outstrength.trait.corrected)
  tmp3 = ANTs:::stat.p(tmp$Outstrength.trait.rand.corrected)
  tmp4 = as.data.frame(rbind(tmp2, tmp3))
  tmp4$type = c("non random", "random")
  tmp4$p = c(tmp$Outstrength.parametric.p.non.rand.corrected[1], tmp$Outstrength.parametric.p.rand.corrected[1])
  result10 = rbind(result10,tmp4)
  
  tmp4 = ANTs:::stat.p(tmp$Degree.trait.corrected)
  tmp5 = ANTs:::stat.p(tmp$Degree.trait.rand.corrected)
  tmp6 = as.data.frame(rbind(tmp4, tmp5))
  tmp6$type = c("non random", "random")
  tmp6$p = c(tmp$Degree.parametric.p.non.rand.corrected[1], tmp$Degree.parametric.p.rand.corrected[1])
  result11 = rbind(result11,tmp6)
  
  tmp6 = ANTs:::stat.p(tmp$Eigen.trait.corrected)
  tmp7 = ANTs:::stat.p(tmp$Eigen.trait.rand.corrected)
  tmp8 = as.data.frame(rbind(tmp6, tmp7))
  tmp8$type = c("non random", "random")
  tmp8$p = c(tmp$Eigen.parametric.p.non.rand.corrected[1], tmp$Eigen.parametric.p.rand.corrected[1])
  result12 = rbind(result12,tmp8)
  
  
  
  cat("#################################################################################", '\n')
  cat("Without GI", '\n')
  cat("Rates of false positives for outstrength : ",
      (nrow(result7[result7$type %in% "random" & result7$`p-value_left_side` < 0.05, ])*100)/nrow(result7[result7$type %in% "random",]), "\n")
  cat("Rates of false negatives for outstrength : ",
      (nrow(result7[result7$type %in% "non random" & result7$`p-value_left_side` > 0.05, ])*100)/nrow(result7[result7$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for outdegree : ",
      (nrow(result8[result8$type %in% "random" & result8$`p-value_left_side` < 0.05, ])*100)/nrow(result8[result8$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree : ",
      (nrow(result8[result8$type %in% "non random" & result8$`p-value_left_side` > 0.05, ])*100)/nrow(result8[result8$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for eigenvector : ",
      (nrow(result9[result9$type %in% "random" & result9$`p-value_left_side` < 0.05, ])*100)/nrow(result9[result9$type %in% "random",]), "\n")
  cat("Rates of false negatives for eigenvector : ",
      (nrow(result9[result9$type %in% "non random" & result9$`p-value_left_side` > 0.05, ])*100)/nrow(result9[result9$type %in% "non random",]), "\n")
  
  cat("\n")
  # Rates of false negatives
  cat("Parametric rates  of false positives for outstrength : ",
      (nrow(result7[result7$type %in% "random" & result7$p < 0.05, ])*100)/nrow(result7[result7$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for outstrength : ",
      (nrow(result7[result7$type %in% "non random" & result7$p > 0.05, ])*100)/nrow(result7[result7$type %in% "non random",]), "\n")
  
  cat("Parametric rates  of false positives for outdegree : ",
      (nrow(result8[result8$type %in% "random" & result8$p < 0.05, ])*100)/nrow(result8[result8$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for outdegree : ",
      (nrow(result8[result8$type %in% "non random" & result8$p > 0.05, ])*100)/nrow(result8[result8$type %in% "non random",]), "\n")
  
  cat("Parametric rates  of false positives for eigenvector : ",
      (nrow(result9[result9$type %in% "random" & result9$p < 0.05, ])*100)/nrow(result9[result9$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for eigenvector : ",
      (nrow(result9[result9$type %in% "non random" & result9$p > 0.05, ])*100)/nrow(result9[result9$type %in% "non random",]), "\n")
  cat("\n")
  
  cat("#################################################################################", '\n')
  cat("With GI", '\n')
  cat("Rates of false positives for outstrength: ",
      (nrow(result10[result10$type %in% "random" & result10$`p-value_left_side` < 0.05, ])*100)/nrow(result10[result10$type %in% "random",]), "\n")
  cat("Rates of false negatives for outstrength: ",
      (nrow(result10[result10$type %in% "non random" & result10$`p-value_left_side` > 0.05, ])*100)/nrow(result7[result10$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for outdegree: ",
      (nrow(result11[result11$type %in% "random" & result11$`p-value_left_side` < 0.05, ])*100)/nrow(result11[result11$type %in% "random",]), "\n")
  cat("Rates of false negatives for outdegree: ",
      (nrow(result11[result11$type %in% "non random" & result11$`p-value_left_side` > 0.05, ])*100)/nrow(result11[result11$type %in% "non random",]), "\n")
  
  cat("Rates of false positives for eigenvector: ",
      (nrow(result12[result12$type %in% "random" & result12$`p-value_left_side` < 0.05, ])*100)/nrow(result12[result12$type %in% "random",]), "\n")
  cat("Rates of false negatives for eigenvector: ",
      (nrow(result12[result12$type %in% "non random" & result12$`p-value_left_side` > 0.05, ])*100)/nrow(result12[result12$type %in% "non random",]), "\n")
  
  cat("\n")
  # Rates of false negatives
  cat("Parametric rates  of false positives for outstrength with GI : ",
      (nrow(result10[result10$type %in% "random" & result10$p < 0.05, ])*100)/nrow(result10[result10$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for outstrength with GI : ",
      (nrow(result10[result10$type %in% "non random" & result10$p > 0.05, ])*100)/nrow(result10[result10$type %in% "non random",]), "\n")
  
  cat("Parametric rates  of false positives for outdegree with GI : ",
      (nrow(result11[result11$type %in% "random" & result11$p < 0.05, ])*100)/nrow(result11[result11$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for outdegree with GI : ",
      (nrow(result11[result11$type %in% "non random" & result11$p > 0.05, ])*100)/nrow(result11[result11$type %in% "non random",]), "\n")
  
  cat("Parametric rates  of false positives for eigenvector with GI : ",
      (nrow(result12[result12$type %in% "random" & result12$p < 0.05, ])*100)/nrow(result12[result12$type %in% "random",]), "\n")
  cat("Parametric rates  of false negatives for eigenvector with GI : ",
      (nrow(result12[result12$type %in% "non random" & result12$p > 0.05, ])*100)/nrow(result12[result12$type %in% "non random",]), "\n")
  cat("\n")
}


#######################################
##### Results
#######################################

d1 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result1[result1$type %in% "random" & result1$`p-value_left_side` < 0.05, ])*100)/nrow(result1[result1$type %in% "random",]),
                    (nrow(result1[result1$type %in% "random" & result1$p < 0.05, ])*100)/nrow(result1[result1$type %in% "random",])),
  "outdegree" = c((nrow(result2[result2$type %in% "random" & result2$`p-value_left_side` < 0.05, ])*100)/nrow(result2[result2$type %in% "random",]),
                  (nrow(result2[result2$type %in% "random" & result2$p < 0.05, ])*100)/nrow(result2[result2$type %in% "random",])),
  "outeigenvector" = c((nrow(result3[result3$type %in% "random" & result3$`p-value_left_side` < 0.05, ])*100)/nrow(result3[result3$type %in% "random",]),
                       (nrow(result3[result3$type %in% "random" & result3$p < 0.05, ])*100)/nrow(result3[result3$type %in% "random",])),
  "Error Type" = rep("False positives rates", 2),
  "Biases" = rep(TRUE, 2),
  "GI" = rep(FALSE, 2)
)

d2 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result1[result1$type %in% "non random" & result1$`p-value_left_side` > 0.05, ])*100)/nrow(result1[result1$type %in% "non random",]),
                    (nrow(result1[result1$type %in% "non random" & result1$p > 0.05, ])*100)/nrow(result1[result1$type %in% "non random",])),
  "outdegree" = c((nrow(result2[result2$type %in% "non random" & result2$`p-value_left_side` > 0.05, ])*100)/nrow(result2[result2$type %in% "non random",]),
                  (nrow(result2[result2$type %in% "non random" & result2$p > 0.05, ])*100)/nrow(result2[result2$type %in% "non random",])),
  "outeigenvector" = c((nrow(result3[result3$type %in% "non random" & result3$`p-value_left_side` > 0.05, ])*100)/nrow(result3[result3$type %in% "non random",]),
                       (nrow(result3[result3$type %in% "non random" & result3$p > 0.05, ])*100)/nrow(result3[result3$type %in% "non random",])),
  "Error Type" = rep("False negatives rates", 2),
  "Biases" = rep(TRUE, 2),
  "GI" = rep(FALSE, 2)
)

d3 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result4[result4$type %in% "random" & result4$`p-value_left_side` < 0.05, ])*100)/nrow(result4[result4$type %in% "random",]),
                    (nrow(result4[result4$type %in% "random" & result4$p < 0.05, ])*100)/nrow(result4[result4$type %in% "random",])),
  "outdegree" = c((nrow(result5[result5$type %in% "random" & result5$`p-value_left_side` < 0.05, ])*100)/nrow(result5[result5$type %in% "random",]),
                  (nrow(result5[result5$type %in% "random" & result5$p < 0.05, ])*100)/nrow(result5[result5$type %in% "random",])),
  "outeigenvector" = c((nrow(result6[result6$type %in% "random" & result6$`p-value_left_side` < 0.05, ])*100)/nrow(result6[result6$type %in% "random",]),
                       (nrow(result6[result6$type %in% "random" & result6$p < 0.05, ])*100)/nrow(result6[result6$type %in% "random",])),
  "Error Type" = rep("False negatives rates", 2),
  "Biases" = rep(TRUE, 2),
  "GI" = rep(TRUE, 2)
)

d4 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result4[result4$type %in% "non random" & result4$`p-value_left_side` > 0.05, ])*100)/nrow(result4[result4$type %in% "non random",]),
                    (nrow(result4[result4$type %in% "non random" & result4$p > 0.05, ])*100)/nrow(result4[result4$type %in% "non random",])),
  "outdegree" = c((nrow(result5[result5$type %in% "non random" & result5$`p-value_left_side` > 0.05, ])*100)/nrow(result5[result5$type %in% "non random",]),
                  (nrow(result5[result5$type %in% "non random" & result5$p > 0.05, ])*100)/nrow(result5[result5$type %in% "non random",])),
  "outeigenvector" = c((nrow(result6[result6$type %in% "non random" & result6$`p-value_left_side` > 0.05, ])*100)/nrow(result6[result6$type %in% "non random",]),
                       (nrow(result6[result6$type %in% "non random" & result6$p > 0.05, ])*100)/nrow(result6[result6$type %in% "non random",])),
  "Error Type" = rep("False negatives rates", 2),
  "Biases" = rep(TRUE, 2),
  "GI" = rep(TRUE, 2)
)


d5 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result7[result7$type %in% "random" & result7$`p-value_left_side` < 0.05, ])*100)/nrow(result7[result7$type %in% "random",]),
                    (nrow(result7[result7$type %in% "random" & result7$p < 0.05, ])*100)/nrow(result7[result7$type %in% "random",])),
  "outdegree" = c((nrow(result8[result8$type %in% "random" & result8$`p-value_left_side` < 0.05, ])*100)/nrow(result8[result8$type %in% "random",]),
                  (nrow(result8[result8$type %in% "random" & result8$p < 0.05, ])*100)/nrow(result8[result8$type %in% "random",])),
  "outeigenvector" = c((nrow(result9[result9$type %in% "random" & result9$`p-value_left_side` < 0.05, ])*100)/nrow(result9[result9$type %in% "random",]),
                       (nrow(result9[result9$type %in% "random" & result9$p < 0.05, ])*100)/nrow(result9[result9$type %in% "random",])),
  "Error Type" = rep("False positives rates", 2),
  "Biases" = rep(FALSE, 2),
  "GI" = rep(FALSE, 2)
)

d6 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result7[result7$type %in% "non random" & result7$`p-value_left_side` > 0.05, ])*100)/nrow(result7[result7$type %in% "non random",]),
                    (nrow(result7[result7$type %in% "non random" & result7$p > 0.05, ])*100)/nrow(result7[result7$type %in% "non random",])),
  "outdegree" = c((nrow(result8[result8$type %in% "non random" & result8$`p-value_left_side` > 0.05, ])*100)/nrow(result8[result8$type %in% "non random",]),
                  (nrow(result8[result8$type %in% "non random" & result8$p > 0.05, ])*100)/nrow(result8[result8$type %in% "non random",])),
  "outeigenvector" = c((nrow(result9[result9$type %in% "non random" & result9$`p-value_left_side` > 0.05, ])*100)/nrow(result9[result9$type %in% "non random",]),
                       (nrow(result9[result9$type %in% "non random" & result9$p > 0.05, ])*100)/nrow(result9[result9$type %in% "non random",])),
  "Error Type" = rep("False negatives rates", 2),
  "Biases" = rep(FALSE, 2),
  "GI" = rep(FALSE, 2)
)

d7 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result10[result10$type %in% "random" & result10$`p-value_left_side` < 0.05, ])*100)/nrow(result10[result10$type %in% "random",]),
                    (nrow(result10[result10$type %in% "random" & result10$p < 0.05, ])*100)/nrow(result10[result10$type %in% "random",])),
  "outdegree" = c((nrow(result11[result11$type %in% "random" & result11$`p-value_left_side` < 0.05, ])*100)/nrow(result11[result11$type %in% "random",]),
                  (nrow(result11[result11$type %in% "random" & result11$p < 0.05, ])*100)/nrow(result11[result11$type %in% "random",])),
  "outeigenvector" = c((nrow(result12[result12$type %in% "random" & result12$`p-value_left_side` < 0.05, ])*100)/nrow(result12[result12$type %in% "random",]),
                       (nrow(result12[result12$type %in% "random" & result12$p < 0.05, ])*100)/nrow(result12[result12$type %in% "random",])),
  "Error Type" = rep("False negatives rates", 2),
  "Biases" = rep(FALSE, 2),
  "GI" = rep(TRUE, 2)
)

d8 = data.frame(
  "approches" = c("Nertwork permutations", "Parametric"),
  "outstrength" = c((nrow(result10[result10$type %in% "non random" & result10$`p-value_left_side` > 0.05, ])*100)/nrow(result10[result10$type %in% "non random",]),
                    (nrow(result10[result10$type %in% "non random" & result10$p > 0.05, ])*100)/nrow(result10[result10$type %in% "non random",])),
  "outdegree" = c((nrow(result11[result11$type %in% "non random" & result11$`p-value_left_side` > 0.05, ])*100)/nrow(result11[result11$type %in% "non random",]),
                  (nrow(result11[result11$type %in% "non random" & result11$p > 0.05, ])*100)/nrow(result11[result11$type %in% "non random",])),
  "outeigenvector" = c((nrow(result12[result12$type %in% "non random" & result12$`p-value_left_side` > 0.05, ])*100)/nrow(result12[result12$type %in% "non random",]),
                       (nrow(result12[result12$type %in% "non random" & result12$p > 0.05, ])*100)/nrow(result12[result12$type %in% "non random",])),
  "Error Type" = rep("False negatives rates", 2),
  "Biases" = rep(FALSE, 2),
  "GI" = rep(TRUE, 2)
)


RESULTS = rbind(d1, d2, d3, d4, d5, d6, d7, d8)

write.csv(RESULTS, file = "results simulation 4.csv")
