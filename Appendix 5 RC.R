#This code is the original code from Farine & G 2020 (https://www.biorxiv.org/content/biorxiv/early/2020/08/04/2020.08.02.232710.full.pdf) to simalte a space confounding variable between individuals sociality and a continuous trait.
#Code have bee warp within a function to call it eseally for each scenarios and each social metrics.
# Simulation of space confounding variable -------------
#library(asnipe)
library(sna)
library(ANTs)
library(ggplot2)
library(ggpubr)
sim.data <- function(N =50, n.obs = 100, effect = T, metric = "DEGREE"){
  
  # create a data frame with IDs
  ids <- data.frame(ID=1:N)
  
  # Generate some possible groups
  groups <- 0.5*N*n.obs
  
  # individual trait values
  ids$TRAIT <- rnorm(N,0,2)
  
  # make these all positive
  ids$TRAIT <- ids$TRAIT-min(ids$TRAIT) + 1
  
  
  # assign number of observations per individual (sums to N*n.obs)
  ids$OBS <- rpois(N,n.obs)
  m <- N*n.obs
  d <- sum(ids$OBS)-m
  if (d > 0) {
    rem <- sample(1:N,d,prob=ids$OBS,replace=TRUE)
    rem <- table(rem)
    ids$OBS[as.numeric(names(rem))] <- ids$OBS[as.numeric(names(rem))]-rem
  }
  if (d < 0) {
    add <- sample(1:N,abs(d),prob=ids$OBS,replace=TRUE)
    add <- table(add)
    ids$OBS[as.numeric(names(add))] <- ids$OBS[as.numeric(names(add))]+add
  }
  ids$OBS[ids$OBS <= 0] <- 1
  
  # create a vector of group sizes (larger number = more individuals)
  group_size <- sample(c(1:10),groups,replace=TRUE)
  
  # Create blank observation matrix
  obs <- matrix(0,nrow=groups,ncol=N)
  
  
  # Allocate individuals to groups, giving smaller trait values a higher probability to be in smaller groups
  if (effect) {
    which.ind <- order(ids$TRAIT, decreasing=FALSE) # allocate in increasing order of TRAIT
  } else {
    which.ind <- 1:nrow(ids) 	# if random, simply select individuals in random order
  }
  
  for (i in which.ind) {
    
    group_size.tmp <- group_size-rowSums(obs)
    probs <- 1/(1+(group_size.tmp)^2)
    probs[group_size.tmp==0] <- 0
    g <- sample(1:groups,ids$OBS[i],prob=probs)
    obs[g,i] <- 1
    
  }
  
  # remove empty groups
  if (any(rowSums(obs) == 0)) {
    group_size <- group_size[-which(rowSums(obs)==0)]
    obs <- obs[-which(rowSums(obs)==0),]
  }
  
  # calculate most observed location per ind
  max_loc <- function(groups,group_size) {
    grps <- table(group_size[groups==1])
    max <- names(grps)[which(grps==max(grps))]
    if(length(max) > 1){
      max <- sample(max, 1)
    }
    return(max)
  }
  
  ids$LOCATION <- apply(obs,2,max_loc,group_size)
  ids$LOCATION.strength =  ids$LOCATION 
  ids$LOCATION.eigen =  ids$LOCATION 
  ids$LOCATION.degree =  ids$LOCATION 
  ids$OBS <- colSums(obs)
  ids$OBS.strength = ids$OBS
  ids$OBS.eigen  = ids$OBS
  ids$OBS.degree = ids$OBS
  
  
  ## NOW DO NETWORK ANALYSIS ON THESE DATA
  
  # Calculate network
  # Removing asnipe for faster computations speed
  #network <- get_network(obs)
  network <- assoc.indices(obs)
  network.corrected <- network/ANTs:::time.heterogeneity(ids$OBS)
  network.corrected[is.nan(network.corrected)] = 0
  network.corrected[is.infinite(network.corrected)] = 0
  ######################################################################################################
  ############ Modification 1
  ######################################################################################################

  
  #get_network_metric <- function(network, metric) {# Calculate Metric
  #  if (metric == "DEGREE") {
  #    METRIC <- rowSums(network)
  #  } else if (metric == "EVCENT") {
  #    METRIC <- sna::evcent(network)
  #  } else if (metric == "BETWEEN") {
  #    METRIC <- sna::betweenness(1/network, ignore.eval=FALSE)
  #  } else {
  #    stop("INCORRECT METRIC SPECIFIED")
  #  }
  #  return(METRIC)
  #}
  #ids$METRIC <- get_network_metric(network, metric)
  ######################################################################################################
  ############ Modification 2
  ######################################################################################################
  ids$strength <- rowSums(network)
  p1 = ggplot(ids, aes(x = LOCATION, y = strength))+geom_point()
  
  ids$strength.corrected <- rowSums(network.corrected)
  p2 = ggplot(ids, aes(x = LOCATION, y = strength.corrected))+geom_point()
  
  ids$eigen <- sna::evcent(network)
  ids$eigen.corrected <- sna::evcent(network.corrected)
  
  ids$degree <- met.degree(network)
  ids$degree.corrected <- met.degree(network.corrected)/ids$OBS
  
  ### Fix effect relationship (mostly for betweenness which is hard to model)
  # do this by sampling individuals in order from largest to smallest network metric
  # and select traits proportionally to their number
  if (effect & cor(ids$TRAIT,ids$strength) < 0.2) {
    metric.order <- order(ids$strength,decreasing=FALSE)
    metric.order <- metric.order[order(rnorm(length(metric.order),mean=1:length(metric.order),sd=length(metric.order)^(1/3)))]
    ids <- ids[order(ids$TRAIT,decreasing=FALSE),]
    obs <- obs[,metric.order]
    network <- assoc.indices(obs)
    ids$LOCATION.strength <- apply(obs,2,max_loc,group_size)
    ids$OBS.strength <- colSums(obs)
    ids$strength <- rowSums(network)
    
    network.corrected <- network/ANTs:::time.heterogeneity(ids$OBS.degree)
    network.corrected[is.nan(network.corrected)] = 0
    network.corrected[is.infinite(network.corrected)] = 0
    ids$strength.corrected <- rowSums(network.corrected)
  }
  
  if (effect & cor(ids$TRAIT,ids$eigen) < 0.2) {
    metric.order <- order(ids$eigen,decreasing=FALSE)
    metric.order <- metric.order[order(rnorm(length(metric.order),mean=1:length(metric.order),sd=length(metric.order)^(1/3)))]
    ids <- ids[order(ids$TRAIT,decreasing=FALSE),]
    obs <- obs[,metric.order]
    network <- assoc.indices(obs)
    ids$LOCATION.eigen <- apply(obs,2,max_loc,group_size)
    ids$OBS.eigen <- colSums(obs)
    ids$eigen <-  sna::evcent(network)
    
    network.corrected <- network/ANTs:::time.heterogeneity(ids$OBS.eigen)
    network.corrected[is.nan(network.corrected)] = 0
    network.corrected[is.infinite(network.corrected)] = 0
    ids$eigen.corrected <- sna::evcent(network.corrected)
  }
  test = cor(ids$TRAIT,ids$degree)
  # Some simulations return fully connected network
  if(is.na(test)){
    metric.order <- order(ids$degree,decreasing=FALSE)
    metric.order <- metric.order[order(rnorm(length(metric.order),mean=1:length(metric.order),sd=length(metric.order)^(1/3)))]
    ids <- ids[order(ids$TRAIT,decreasing=FALSE),]
    obs <- obs[,metric.order]
    network <- assoc.indices(obs)
    ids$LOCATION.degree <- apply(obs,2,max_loc,group_size)
    ids$OBS.degree <- colSums(obs)
    ids$degree <-  met.degree(network)
    
    network.corrected <- network/ANTs:::time.heterogeneity(ids$OBS.degree)
    network.corrected[is.nan(network.corrected)] = 0
    network.corrected[is.infinite(network.corrected)] = 0
    ids$degree.corrected <- met.degree(network.corrected)
  }else{
    if (effect & cor(ids$TRAIT,ids$degree) < 0.2) {
      metric.order <- order(ids$degree,decreasing=FALSE)
      metric.order <- metric.order[order(rnorm(length(metric.order),mean=1:length(metric.order),sd=length(metric.order)^(1/3)))]
      ids <- ids[order(ids$TRAIT,decreasing=FALSE),]
      obs <- obs[,metric.order]
      network <- assoc.indices(obs)
      ids$LOCATION.degree <- apply(obs,2,max_loc,group_size)
      ids$OBS.degree <- colSums(obs)
      ids$degree <-  met.degree(network)
      
      network.corrected <- network/ANTs:::time.heterogeneity(ids$OBS.degree)
      network.corrected[is.nan(network.corrected)] = 0
      network.corrected[is.infinite(network.corrected)] = 0
      ids$degree.corrected <- met.degree(network.corrected)
    }
  }

  return(ids)
}

# Set parameters
N <- c(10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120)
n.obs <- c(5, 10, 15, 20, 25, 30, 35, 40)
params <- expand.grid(N=N,n.obs=n.obs)

## Effect = T and location = T --------
#he simulated data should have a strong spatial ‘nuisance’ effect, 
#with the local density of individuals varying across space, 
# and the proportion of significant results should again approach 5%.
result0.0 = result0.1 = result0.2 = result0.3 = result0.4 = result0.5 = NULL
result1 = result2 = result3 = result4 = result5 = result6 = NULL
effect = T
location = T
a = 1 
# For each params in values run simulation
for (a in a:nrow(params)) {
  # 100 times
  for(b in 1:100){
    cat("Simulation", a, " variant ", b, "\n")
    dat = sim.data(N = params[a,1], n.obs = params[a,2],effect = effect)
    ## Observations residuals correction
    m1 = lm(strength.corrected ~ LOCATION.strength, data = dat)
    dat$new.strength = resid(m1)
    
    m1 = lm(eigen.corrected ~ LOCATION.eigen, data = dat)
    dat$new.eigen = resid(m1)
    
    m1 = lm(degree.corrected ~ LOCATION.degree, data = dat)
    dat$new.degree = resid(m1)
    
    
    
    test0 = test0.1 = test1 = test2 = test3 = test4 = test5 = test6 = NULL
    for(d in 1:1000){
      if(d==1){
        test0.0 = summary(lm(data = dat, formula = new.strength ~ TRAIT))$coefficients[2,4]
        test0.1 = summary(lm(data = dat, formula = new.eigen ~ TRAIT))$coefficients[2,4]
        test0.2 = summary(lm(data = dat, formula = new.degree ~ TRAIT))$coefficients[2,4]
        
        test0.3 = summary(lm(data = dat, formula = strength.corrected ~ TRAIT))$coefficients[2,4]
        test0.4 = summary(lm(data = dat, formula = eigen.corrected ~ TRAIT))$coefficients[2,4]
        test0.5 = summary(lm(data = dat, formula = degree.corrected ~ TRAIT))$coefficients[2,4]
        
        test1 = rbind(test1, coefficients(lm(data = dat, formula = new.strength ~ TRAIT)))
        test2 = rbind(test2, coefficients(lm(data = dat, formula = new.eigen ~ TRAIT)))
        test3 = rbind(test3, coefficients(lm(data = dat, formula = new.degree ~ TRAIT)))
        
        test4 = rbind(test4, coefficients(lm(data = dat, formula = strength.corrected ~ TRAIT)))
        test5 = rbind(test5, coefficients(lm(data = dat, formula = eigen.corrected ~ TRAIT)))
        test6 = rbind(test6, coefficients(lm(data = dat, formula = degree.corrected ~ TRAIT)))

        }else{
        test1 = rbind(test1,coefficients(lm(data = dat, formula = new.strength ~ sample(TRAIT))))
        test2 = rbind(test2,coefficients(lm(data = dat, formula = new.eigen ~ sample(TRAIT))))
        test3 = rbind(test3,coefficients(lm(data = dat, formula = new.degree ~ sample(TRAIT))))
        
        test4 = rbind(test4,coefficients(lm(data = dat, formula = strength.corrected ~ sample(TRAIT))))
        test5 = rbind(test5,coefficients(lm(data = dat, formula = eigen.corrected ~ sample(TRAIT))))
        test6 = rbind(test6,coefficients(lm(data = dat, formula = degree.corrected ~ sample(TRAIT))))
        
      }
    }
    
    #warp information
    result0.0 =  rbind(result0.0,data.frame("metric" = "strength", "p" = test0.0,N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result0.1 =  rbind(result0.1,data.frame("metric" = "eigen", "p" =test0.1,N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result0.2 =  rbind(result0.2,data.frame("metric" = "degree", "p" =test0.2,N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    
    result0.3 =  rbind(result0.3,data.frame("metric" = "strength", "p" =test0.3,N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result0.4 =  rbind(result0.4,data.frame("metric" = "eigen", "p" =test0.4,N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result0.5 =  rbind(result0.5,data.frame("metric" = "degree", "p" =test0.5,N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    
    result1 =  rbind(result1,data.frame("metric" = "strength", t(ANTs:::stat.p(test1[,2])),N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result2 =  rbind(result2,data.frame("metric" = "eigen", t(ANTs:::stat.p(test2[,2])),N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result3 =  rbind(result3,data.frame("metric" = "degree", t(ANTs:::stat.p(test3[,2])),N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    
    result4 =  rbind(result4,data.frame("metric" = "strength", t(ANTs:::stat.p(test4[,2])),N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result5 =  rbind(result5,data.frame("metric" = "eigen", t(ANTs:::stat.p(test5[,2])),N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    result6 =  rbind(result6,data.frame("metric" = "degree", t(ANTs:::stat.p(test6[,2])),N = params[a,1], n.obs = params[a,2], "effect" = effect, "location" = location))
    
    cat("Parametric false positives rates with residual correction for strength: ", (nrow(result0.0[result0.0$p < 0.05,])*100)/nrow(result0.0), "\n")
    cat("Parametric false positives rates with residual correction for eigen: ", (nrow(result0.1[result0.1$p < 0.05,])*100)/nrow(result0.1), "\n")
    cat("Parametric false positives rates with residual correction for degrtee: ", (nrow(result0.2[result0.2$p < 0.05,])*100)/nrow(result0.2), "\n")
    
    cat("Parametric false positives rates without residual correction for strength: ", (nrow(result0.3[result0.3$p < 0.05,])*100)/nrow(result0.3), "\n")
    cat("Parametric false positives rates without residual correction for eigen: ", (nrow(result0.4[result0.4$p < 0.05,])*100)/nrow(result0.4), "\n")
    cat("Parametric false positives rates without residual correction for degrtee: ", (nrow(result0.5[result0.5$p < 0.05,])*100)/nrow(result0.5), "\n")    
    
    cat("\n")
    cat("False positives rates with residual correction for strength: ", (nrow(result1[result1$p.value_left_side < 0.05,])*100)/nrow(result1), "\n")
    cat("False positives rates with residual correction for eigen: ", (nrow(result2[result2$p.value_left_side < 0.05,])*100)/nrow(result2), "\n")
    cat("False positives rates with residual correction for degrtee: ", (nrow(result3[result3$p.value_left_side < 0.05,])*100)/nrow(result3), "\n")

    cat("False positives rates without residual correction for strength: ", (nrow(result4[result4$p.value_left_side < 0.05,])*100)/nrow(result4), "\n")
    cat("False positives rates without residual correction for eigen: ", (nrow(result5[result5$p.value_left_side < 0.05,])*100)/nrow(result5), "\n")
    cat("False positives rates without residual correction for degrtee: ", (nrow(result6[result6$p.value_left_side < 0.05,])*100)/nrow(result6), "\n")    
  }
}

