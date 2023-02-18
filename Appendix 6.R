#' @param N an integer indicating the total number of individuals
#' @param Obs.mean.bias amount of observation bias between sex (in favor of females)
#' @param diff.soc.pheno amount of interaction differences per observation between sex (in favor of females)
#' @param N.obs.pheno the base number of observation per individuals (will be define for males)
#' @param N.alters.per.obs the number of alters per observation (will be define for males)
sim <- function(N, Obs.mean.bias, diff.soc.pheno, N.obs.pheno, N.alters.per.obs){
  ids = 1 : N
  sex = sample(c("F", "M"), N, replace = T, prob = c(0.5,0.5))

  # observations biases--------------
  obs.m = N.obs.pheno
  obs.f = obs.m + ((obs.m*Obs.mean.bias)/100) # add the amount of biases define by "Obs.mean.bias"

  # soc biases
  N.m = N.alters.per.obs
  N.f = N.alters.per.obs + ((N.alters.per.obs*diff.soc.pheno)/100)# add the amount of interactions per observation define by "diff.soc.pheno"

  Result = NULL
  sampling.effort = NULL
  for (a in 1:N) {
    ego = a
    ego.sex = sex[a]
    if(ego.sex == "M"){
      ego.obs = rnorm(1, mean = obs.m, sd = (obs.m*1)/100)# 1% of observation variation within phenotype
      ego.obs = round(ego.obs)
      sampling.effort = rbind(sampling.effort, data.frame(ego, "sampling.effort" = ego.obs))

      for (b in 1:ego.obs) {
        cat("Processing individual ", a, " observation #", b, "\r")
        alters = sample(ids[!ids %in% ego], rnorm(1, mean = N.m, sd = (N.m*5)/100), replace = T)
        Result = rbind(Result, data.frame(
          ego, "sex" = ego.sex, "focal" = b, "alters" = alters
        ))
      }

    }else{
      ego.obs = rnorm(1, mean = obs.f, sd = (obs.f*1)/100)# 10% of observation variation within phenotype
      ego.obs = round(ego.obs)
      sampling.effort = rbind(sampling.effort, data.frame(ego, "sampling.effort" = ego.obs))

      for (b in 1:ego.obs) {
        cat("Processing individual ", a, " observation #", b, "\r")
        alters = sample(ids[!ids %in% ego], round(rnorm(1, mean = N.f, sd = (N.f*10)/100)), replace = T)
        Result = rbind(Result, data.frame(
          ego, "sex" = ego.sex, "focal" = b, "alters" = alters
        ))
      }
    }
  }
  return(list(Result, ids, sex, sampling.effort))
}
sim(N = 30, Obs.mean.bias = 0, diff.soc.pheno = 0, N.obs.pheno= 50, N.alters.per.obs = 10, print = T)

library(ANTs)
library(ggplot2)
library(ggpubr)
library(sjPlot)
SIM <- function(N = 30, Obs.mean.bias = -50, diff.soc.pheno = -10, N.obs.pheno = 50, N.alters.per.obs = 10, print = F){
  Result = sim(N = N, Obs.mean.bias =Obs.mean.bias, diff.soc.pheno = diff.soc.pheno, N.obs.pheno= N.obs.pheno, N.alters.per.obs = N.alters.per.obs)
  m = df.to.mat(Result[[1]], actor = "ego", receiver = "alters", num.ids = T)
  data = data.frame("ids" = Result[[2]], "sex" = Result[[3]])
  data$sampling.effort = Result[[4]]$sampling.effort
  p0 = ggplot(data, aes(x = sex, y = sampling.effort, group = sex))+geom_boxplot()+geom_point()

  # Ivan correction----------------------------------------
  m2 = m / data$sampling.effort
  m2 = (m2 + t(m2))/2

  data$strength.NEW = met.strength(m2)
  p1 = ggplot(data, aes(x = sex, y = strength.NEW, group = sex))+geom_boxplot()+geom_point()


  # No correction----------------------------------------
  m = df.to.mat(Result[[1]], actor = "ego", receiver = "alters", num.ids = T)
  #data = data.frame("ids" = Result[[2]], "sex" = Result[[3]])
  #data$sampling.effort = Result[[4]]$sampling.effort
  data$strength.no.cor = met.strength(m)
  p1 = ggplot(data, aes(x = sex, y = strength.no.cor, group = sex))+geom_boxplot()+geom_point()

  # Correction by frequencies ----------------------------------------
  m = df.to.mat(Result[[1]], actor = "ego", receiver = "alters", num.ids = T, tobs = Result[[4]]$sampling.effort)
  data$strength.cor.fr = met.strength(m)
  p2 = ggplot(data, aes(x = sex, y = strength.cor.fr, group = sex))+geom_boxplot()+geom_point()


  # Correction by indices----------------------------------------
  m = df.to.mat(Result[[1]], actor = "ego", receiver = "alters", num.ids = T)
  fa = matrix(0, ncol = nrow(data), nrow = nrow(data))
  for (a in 1:nrow(data)) {
    fa[a,] = data$sampling.effort[a]
  }

  fb = matrix(0, ncol = nrow(data), nrow = nrow(data))
  for (a in 1:nrow(data)) {
    fb[a,] = data$sampling.effort
  }

  ya = abs(fa - m)

  yb = abs(fb - m)


  sri <- ((m) /(m + ya + yb ))
  diag(sri) = 0
  data$strength.cor.assoc = met.strength(sri)

  p3 = ggplot(data, aes(x = sex, y = strength.cor.assoc, group = sex))+geom_boxplot()+geom_point()

  # Correction by GI----------------------------------------
  tmp = Result[[1]]
  obs = unlist(lapply(split(tmp, tmp$ego), nrow))
  obs = obs/data$sampling.effort
  mean.obs = mean(obs, na.rm = T)
  dif.mean = obs - mean.obs
  dif.mean = abs(dif.mean)
  dif.mean = ifelse(dif.mean == 0 , 0.000001, dif.mean)

  data$dif.mean = obs
  ggplot(data, aes(x = sex, y = dif.mean))+ geom_boxplot()+geom_point()

  fa = matrix(0, ncol = nrow(data), nrow = nrow(data))
  for (a in 1:nrow(data)) {
    fa[a,] = dif.mean[a]
  }

  fb = t(fa)


  sampling.effort = (fa + fb)
  GI = (sri)/(sampling.effort)

  diag(GI)=0
  data$strength.cor.GI = met.strength(GI)




  tmp = Result[[1]]
  tmp$unique = paste(tmp$ego, tmp$focal)
  tmp = split(tmp,tmp$unique)
  t = lapply(tmp, function(x){
    c(unique(x$ego), unique(x$sex), nrow(x))
  })
  t = do.call("rbind", t)
  t = as.data.frame(t)
  t$V3 = as.numeric(t$V3)
  colnames(t) = c("ego", "sex", "alters")

  ggplot(t, aes( x = sex, y = alters))+geom_boxplot()+geom_point()

  unlist(lapply(split(t, t$ego), nrow))
  #plot(as.vector(GI),as.vector(sri))
  #plot(as.vector(sampling.effort),as.vector(sri))
  #plot(as.vector(sampling.effort),as.vector(GI))

  #tmp = data.frame(as.vector(GI), as.vector(sri), as.vector(sampling.effort))
  #colnames(tmp) = c("gi", "sri", "sampling")
  #ggplot(tmp, aes(x = gi, y = sri, color = sampling))+geom_point()
  sri[sri == 0] = 0.000000001

  #plot(log(GI))
  #plot((GI))
  #plot(log(GI^(1/3)))
  #GI2 = -((sri)/(log(GI))) #Work but not for eigenvector and 0 can't be handle
  #GI = (1/sri)/GI
  #GI = sri/(1/GI)
  #GI = GI
  diag(GI) = 0

  GI[is.infinite(GI)] = 0
  GI[is.nan(GI)] = 0
  data$strength.cor.GI = met.strength(GI)

  p4 = ggplot(data, aes(x = sex, y = strength.cor.GI, group = sex))+geom_boxplot()+geom_point()

  # regression-----
  data$strength.cor.resid =resid(lm(data = data, formula = strength.no.cor ~ sampling.effort))
  p5 = ggplot(data, aes(x = sex, y = strength.cor.resid, group = sex))+geom_boxplot()+geom_point()

  #Results------------
  if(print){
    print(ggarrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3))
  }

  m1 = lm(data = data, formula = strength.no.cor ~ sex)
  m2 = lm(data = data, formula = strength.cor.fr ~ sex)
  m3 = lm(data = data, formula = strength.cor.assoc ~ sex)
  m4 = lm(data = data, formula = strength.cor.GI ~ sex)
  m5 = lm(data = data, formula = strength.cor.resid ~ sex)
  tab_model(m1,m2,m3,m4,m5)
  t = c(summary(m1)$coefficients[2,4],
           summary(m2)$coefficients[2,4],
           summary(m3)$coefficients[2,4],
           summary(m4)$coefficients[2,4],
           summary(m5)$coefficients[2,4])
  return(data)
}
d = SIM(Obs.mean.bias = 0, diff.soc.pheno = 0, print = T)
summary(lm(strength.cor.GI~sex, data = d))

diff.soc.pheno <- seq(from = -50, to = 50, by = 10)
Obs.mean.bias <- -50
N <- 50
N.obs.pheno <- 50
N.alters.per.obs = 10
params <- expand.grid(diff.soc.pheno = diff.soc.pheno,Obs.mean.bias = Obs.mean.bias, N = N, N.obs.pheno = N.obs.pheno, N.alters.per.obs = N.alters.per.obs )
nrow(params)
RESULT = NULL
a = 1
for (a in a:nrow(params)) {
  cat("#################################################################################", '\n')
  cat("SIM ", a, "/", nrow(params), ", with Obs.mean.bias = ", params$Obs.mean.bias[a], ", diff.soc.pheno = ", params$diff.soc.pheno[a], ", N = ", params$N[a],
      ", N.obs.pheno = ",  params$N.obs.pheno[a], ", N.alters.per.obs = ", params$N.alters.per.obs[a], "\n")
  tmp =SIM(N = params$N[a], Obs.mean.bias = params$Obs.mean.bias[a], diff.soc.pheno = params$diff.soc.pheno[a],
                            N.obs.pheno= params$N.obs.pheno[a], N.alters.per.obs = params$N.alters.per.obs[a])
  tmp$sim = paste(params$diff.soc.pheno[a], sep ="","/", params$Obs.mean.bias[a])
  tmp$sim = as.factor(tmp$sim)
  tmp$diff.soc.pheno =  params$diff.soc.pheno[a]
  tmp$Obs.mean.bias = params$Obs.mean.bias[a]
  RESULT = rbind(RESULT, tmp)

  print(ggplot(RESULT, aes(x = sim , y = strength.cor.fr, fill= sex))+geom_boxplot())+facet_grid(~diff.soc.pheno)
}
colnames(RESULT)
print(ggplot(RESULT, aes(x = sim , y = strength.cor.fr, fill= sex))+geom_boxplot())
print(ggplot(RESULT, aes(x = sim , y = strength.cor.assoc, fill= sex))+geom_boxplot())
print(ggplot(RESULT, aes(x = sim , y = strength.cor.GI, fill= sex))+geom_boxplot()).
print(ggplot(RESULT, aes(x = sim , y = strength.cor.GI2, fill= sex))+geom_boxplot())
print(ggplot(RESULT, aes(x = sim , y = strength.NEW, fill= sex))+geom_boxplot())
