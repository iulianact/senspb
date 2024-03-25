
rm(list=ls())

set.seed(1)

library(sandwich)

source("sim_functions.R")

#-----------------------------------------------------------
nrep <- 1000

N <- 500
pX <- 0.5

b0 <- 0.2; bZ <- 1
alpha_fun <- function(z){b0+bZ*z} 
alpha_fun0 <- function(z){0} 

params_scen1 <- c(a0=1,aX=1,aZ=0)
params_scen2 <- c(a0=1,aX=0,aZ=-1)
params_scen3 <- c(a0=1,aX=2,aZ=-2)

deltagrid <- seq(-5,5,0.1)

form <- as.formula(Y~X+Z)
formis <- as.formula(Y~X+log(abs(Z)))

lev <- 0.05

aggres <- function(reslist,pitrue,zalpha2=qnorm(lev/2,lower.tail=FALSE)){
  pis <- numeric(nrep)
  sds <- numeric(nrep)
  betas <- numeric(nrep)
  coverage <- numeric(nrep)
  cilow <- numeric(nrep)
  ciup <- numeric(nrep)
  
  for(i in 1:nrep){
    pis[i] <- reslist[[i]]$piest
    sds[i] <- sqrt(reslist[[i]]$varpi)
    betas[i] <- reslist[[i]]$betaest
    cilow[i] <- pis[i] - zalpha2*sds[i]
    ciup[i] <- pis[i] + zalpha2*sds[i]
    coverage[i] <- (pitrue < ciup[i]) & (pitrue > cilow[i])
  }
  
  mean_pi <- mean(pis)
  sd_pi <- sd(pis)
  mean_est_sd <- mean(sds)
  mean_beta <- mean(betas)
  sd_beta <- sd(betas)
  mean_coverage <- mean(coverage)
  all <- data.frame(pis=pis,sds=sds,cilow=cilow,ciup=ciup,coverage=coverage,betas=betas)
  
  list(all = all,
       mean_pi=mean_pi,
       sd_pi=sd_pi,
       mean_est_sd=mean_est_sd,
       mean_coverage=mean_coverage,
       mean_beta=mean_beta,
       sd_beta=sd_beta)
}

aggrscen <- function(res10,res01,res11,res00,pitrue,zalpha2=qnorm(lev/2,lower.tail=FALSE)){
  resdf <- data.frame(matrix(NA,nrow=4,ncol=6))
  colnames(resdf) <- c("mean_pi","se_pi","mean_est_se","mean_coverage","mean_beta","sd_beta")
  rownames(resdf) <- c("pi10","pi01","pi11","pi00")
  
  resdf[1,] <- c(res10$mean_pi,res10$sd_pi,res10$mean_est_sd,res10$mean_coverage,res10$mean_beta,res10$sd_beta)
  resdf[2,] <- c(res01$mean_pi,res01$sd_pi,res01$mean_est_sd,res01$mean_coverage,res01$mean_beta,res01$sd_beta)
  resdf[3,] <- c(res11$mean_pi,res11$sd_pi,res11$mean_est_sd,res11$mean_coverage,res11$mean_beta,res11$sd_beta)
  resdf[4,] <- c(res00$mean_pi,res00$sd_pi,res00$mean_est_sd,res00$mean_coverage,res00$mean_beta,res00$sd_beta)
  
  resdf
}


#-----------------------------------------------------------
# Scenario I
#-----------------------------------------------------------
gen_scen1 <- lapply(1:nrep, FUN=function(x){generate(N=N,pX=pX,
                                                     a0=params_scen1["a0"],aX=params_scen1["aX"],aZ=params_scen1["aZ"],
                                                     alpha_fun=alpha_fun0,
                                                     muZ=0,sigZ=1)})

pitrue_scen1 <- gen_scen1[[1]]$p01
betatrue_scen1 <- c(gen_scen1[[1]]$beta10,gen_scen1[[1]]$beta01,gen_scen1[[1]]$beta11,gen_scen1[[1]]$beta00)
deltatrue_scen1 <- c(gen_scen1[[1]]$delta10,gen_scen1[[1]]$delta01,gen_scen1[[1]]$delta11,gen_scen1[[1]]$delta00)

est_scen1_10 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen1[[x]]$dat,x=1,y=0,deltagrid=gen_scen1[[x]]$delta10,formula=form)})
est_scen1_01 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen1[[x]]$dat,x=0,y=1,deltagrid=gen_scen1[[x]]$delta01,formula=form)})
est_scen1_11 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen1[[x]]$dat,x=1,y=1,deltagrid=gen_scen1[[x]]$delta11,formula=form)})
est_scen1_00 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen1[[x]]$dat,x=0,y=0,deltagrid=gen_scen1[[x]]$delta00,formula=form)})

res_scen1_10 <- aggres(est_scen1_10, pitrue = pitrue_scen1)
res_scen1_01 <- aggres(est_scen1_01, pitrue = pitrue_scen1)
res_scen1_11 <- aggres(est_scen1_11, pitrue = pitrue_scen1)
res_scen1_00 <- aggres(est_scen1_00, pitrue = pitrue_scen1)

res_all1 <- cbind(scenario=rep("Scenario I",4),
                  pitrue=rep(pitrue_scen1,4),
                  aggrscen(res_scen1_10,res_scen1_01,res_scen1_11,res_scen1_00,pitrue = pitrue_scen1),
                  beta_true=betatrue_scen1,
                  delta_true=deltatrue_scen1)
res_all1


#-----------------------------------------------------------
# Scenario II
#-----------------------------------------------------------
gen_scen2 <- lapply(1:nrep, FUN=function(x){generate(N=N,pX=pX,
                                                     a0=params_scen2["a0"],aX=params_scen2["aX"],aZ=params_scen2["aZ"],
                                                     alpha_fun=alpha_fun,
                                                     muZ=0,sigZ=1)})

pitrue_scen2 <- gen_scen2[[1]]$p01
betatrue_scen2 <- c(gen_scen2[[1]]$beta10,gen_scen2[[1]]$beta01,gen_scen2[[1]]$beta11,gen_scen2[[1]]$beta00)
deltatrue_scen2 <- c(gen_scen2[[1]]$delta10,gen_scen2[[1]]$delta01,gen_scen2[[1]]$delta11,gen_scen2[[1]]$delta00)

est_scen2_10 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen2[[x]]$dat,x=1,y=0,deltagrid=gen_scen2[[x]]$delta10,formula=form)})
est_scen2_01 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen2[[x]]$dat,x=0,y=1,deltagrid=gen_scen2[[x]]$delta01,formula=form)})
est_scen2_11 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen2[[x]]$dat,x=1,y=1,deltagrid=gen_scen2[[x]]$delta11,formula=form)})
est_scen2_00 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen2[[x]]$dat,x=0,y=0,deltagrid=gen_scen2[[x]]$delta00,formula=form)})

res_scen2_10 <- aggres(est_scen2_10, pitrue = pitrue_scen2)
res_scen2_01 <- aggres(est_scen2_01, pitrue = pitrue_scen2)
res_scen2_11 <- aggres(est_scen2_11, pitrue = pitrue_scen2)
res_scen2_00 <- aggres(est_scen2_00, pitrue = pitrue_scen2) 

res_all2 <- cbind(scenario=rep("Scenario II",4),
                  pitrue=rep(pitrue_scen2,4),
                  aggrscen(res_scen2_10,res_scen2_01,res_scen2_11,res_scen2_00,pitrue = pitrue_scen2),
                  beta_true=betatrue_scen2,
                  delta_true=deltatrue_scen2)
res_all2


#-----------------------------------------------------------
# Scenario III
#-----------------------------------------------------------
gen_scen3 <- lapply(1:nrep, FUN=function(x){generate(N=N,pX=pX,
                                                     a0=params_scen3["a0"],aX=params_scen3["aX"],aZ=params_scen3["aZ"],
                                                     alpha_fun=alpha_fun,
                                                     muZ=0,sigZ=1)})

pitrue_scen3 <- gen_scen3[[1]]$p01
betatrue_scen3 <- c(gen_scen3[[1]]$beta10,gen_scen3[[1]]$beta01,gen_scen3[[1]]$beta11,gen_scen3[[1]]$beta00)
deltatrue_scen3 <- c(gen_scen3[[1]]$delta10,gen_scen3[[1]]$delta01,gen_scen3[[1]]$delta11,gen_scen3[[1]]$delta00)

est_scen3_10 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3[[x]]$dat,x=1,y=0,deltagrid=gen_scen3[[x]]$delta10,formula=form)})
est_scen3_01 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3[[x]]$dat,x=0,y=1,deltagrid=gen_scen3[[x]]$delta01,formula=form)})
est_scen3_11 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3[[x]]$dat,x=1,y=1,deltagrid=gen_scen3[[x]]$delta11,formula=form)})
est_scen3_00 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3[[x]]$dat,x=0,y=0,deltagrid=gen_scen3[[x]]$delta00,formula=form)})

res_scen3_10 <- aggres(est_scen3_10, pitrue = pitrue_scen3)
res_scen3_01 <- aggres(est_scen3_01, pitrue = pitrue_scen3)
res_scen3_11 <- aggres(est_scen3_11, pitrue = pitrue_scen3)
res_scen3_00 <- aggres(est_scen3_00, pitrue = pitrue_scen3)

res_all3 <- cbind(scenario=rep("Scenario III",4),
                  pitrue=rep(pitrue_scen3,4),
                  aggrscen(res_scen3_10,res_scen3_01,res_scen3_11,res_scen3_00,pitrue = pitrue_scen3),
                  beta_true=betatrue_scen3,
                  delta_true=deltatrue_scen3)
res_all3


#-----------------------------------------------------------
# Scenario IIImis
#-----------------------------------------------------------
gen_scen3mis <- gen_scen3

est_scen3mis_10 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3mis[[x]]$dat,x=1,y=0,deltagrid=gen_scen3mis[[x]]$delta10,formula=formis)})
est_scen3mis_01 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3mis[[x]]$dat,x=0,y=1,deltagrid=gen_scen3mis[[x]]$delta01,formula=formis)})
est_scen3mis_11 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3mis[[x]]$dat,x=1,y=1,deltagrid=gen_scen3mis[[x]]$delta11,formula=formis)})
est_scen3mis_00 <- lapply(1:nrep, FUN=function(x){estpb(dat=gen_scen3mis[[x]]$dat,x=0,y=0,deltagrid=gen_scen3mis[[x]]$delta00,formula=formis)})

res_scen3mis_10 <- aggres(est_scen3mis_10, pitrue = pitrue_scen3)
res_scen3mis_01 <- aggres(est_scen3mis_01, pitrue = pitrue_scen3)
res_scen3mis_11 <- aggres(est_scen3mis_11, pitrue = pitrue_scen3)
res_scen3mis_00 <- aggres(est_scen3mis_00, pitrue = pitrue_scen3)

res_all3mis <- cbind(scenario=rep("Scenario IIImis",4),
                     pitrue=rep(pitrue_scen3,4),
                     aggrscen(res_scen3mis_10,res_scen3mis_01,res_scen3mis_11,res_scen3mis_00,pitrue = pitrue_scen3),
                     beta_true=betatrue_scen3,
                     delta_true=deltatrue_scen3)
res_all3mis



