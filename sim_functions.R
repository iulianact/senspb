
generate <- function(N,pX,
                     a0,aX,aZ,
                     alpha_fun,
                     muZ=0,sigZ=1){ 
  Z <- rnorm(N,muZ,sigZ)
  X <- rbinom(N,1,pX)
  Y <- rbinom(N,1,plogis(a0 + aX * X + aZ * Z))
  
  p0Z <- function(z){plogis(a0 + aX * 0 + aZ * z)}
  p1Z <- function(z){plogis(a0 + aX * 1 + aZ * z)}
  
  integrandZ <- function(z,y0,y1){
    alpha <- alpha_fun(z)
    p0 <- p0Z(z)
    p1 <- p1Z(z)
    Delta <- exp(alpha)^2 * (p0-p1)^2 + 2* exp(alpha) * (p0*(1-p0) + p1*(1-p1)) + (p0+p1-1)^2
    denom <- 2*(exp(alpha)-1)
    num_min <- ((p0+p1)*(exp(alpha)-1)+1 - sqrt(Delta))
    sol1 <- ifelse(is.infinite(exp(alpha)) | is.infinite(exp(alpha)^2),min(p0,p1),
                   ifelse(denom==0, p0*p1, num_min/denom))
    picksol <- sol1
    
    p11 <- picksol
    p10 <- p0-p11
    p01 <- p1-p11
    p00 <- 1-p1-p0+p11
    return((y0*y1* p11 + y0*(1-y1)*p10 + (1-y0)*y1*p01 + (1-y0)*(1-y1)*p00)* dnorm(z,muZ,sigZ))
  }
  
  p11 <- integrate(f=Vectorize(integrandZ),lower=-Inf, upper=Inf, y0=1, y1=1)$value
  p10 <- integrate(f=Vectorize(integrandZ),lower=-Inf, upper=Inf, y0=1, y1=0)$value
  p01 <- integrate(f=Vectorize(integrandZ),lower=-Inf, upper=Inf, y0=0, y1=1)$value
  p00 <- integrate(f=Vectorize(integrandZ),lower=-Inf, upper=Inf, y0=0, y1=0)$value
  
  alpha_true <- log((p11*p00)/(p01*p10))

  p1true <- p11+p01
  p0true <- p10+p11
  
  dat <- data.frame(Z=Z,X=X,Y=Y)
  
  params <- list(N=N,pX=pX,
                 a0=a0,aX=aX,aZ=aZ,
                 alpha_fun=alpha_fun,
                 muZ=muZ,sigZ=sigZ)
  
  betaxyreal <- function(x,y){
    integrand <- function(z,yin,xin){
      yin*plogis(a0 + aX * xin + aZ * z) * dnorm(z,muZ,sigZ) + (1-yin)*(1-plogis(a0 + aX * xin + aZ * z)) * dnorm(z,muZ,sigZ)
    }
    pYyX1minx <- integrate(f=integrand, lower=-Inf,upper=Inf,yin=y,xin=1-x)$value
    bigint <- function(z,xin,yin){
      (1/pYyX1minx) * plogis(a0 + aX * xin + aZ * z) * integrand(z,yin,1-xin)
    }
    qlogis(integrate(f=bigint, lower=-Inf,upper=Inf,yin=y,xin=x)$value)
  }
  
  alphaxyreal <- function(x,y){
    integrand <- function(z,yin,xin){
      yin*plogis(a0 + aX * xin + aZ * z) * dnorm(z,muZ,sigZ) + (1-yin)*(1-plogis(a0 + aX * xin + aZ * z)) * dnorm(z,muZ,sigZ)
    }
    if(x==1 & y==1){numerator <- p11}
    if(x==1 & y==0){numerator <- p01}
    if(x==0 & y==1){numerator <- p11}
    if(x==0 & y==0){numerator <- p10}
    qlogis(numerator/integrate(f=integrand,lower=-Inf,upper=Inf,yin=y,xin=1-x)$value)
  }
  
  deltaxyreal <- function(x,y){
    alphaxyreal(x,y) - betaxyreal(x,y)
  }
  
  list(params=params,
       dat=dat,
       p11=p11,p10=p10,p01=p01,p00=p00,
       alpha_true=alpha_true,
       p1true=p1true,
       p0true=p0true,
       beta11=betaxyreal(1,1),beta10=betaxyreal(1,0),beta01=betaxyreal(0,1),beta00=betaxyreal(0,0),
       alpha11=alphaxyreal(1,1),alpha10=alphaxyreal(1,0),alpha01=alphaxyreal(0,1),alpha00=alphaxyreal(0,0),
       delta11=deltaxyreal(1,1),delta10=deltaxyreal(1,0),delta01=deltaxyreal(0,1),delta00=deltaxyreal(0,0))
}


dplogis <- function(x){
  plogis(x)/(1+exp(x))
}


estpb <- function(dat,x,y,deltagrid,formula){
  N <- nrow(dat)
  fit <- glm(formula, data=dat, family="binomial")
  des <- model.matrix(fit)
  res <- residuals(fit, type="response")
  K <- length(coef(fit))
  newdatx <- dat; newdatx$X <- x
  predx <- predict(object=fit,newdata=newdatx,type="response")
  
  q <- mean(dat$X==1-x & dat$Y==y)
  r <- mean(dat$X==1-x)
  betaest <- qlogis(1/nrow(dat)*sum(as.numeric(dat$Y==y & dat$X==1-x)*predx)/q)
  
  if(x==1 & y==0){ piest <- mean(dat$Y[dat$X==0]==0)*plogis(deltagrid+betaest) }
  if(x==0 & y==1){ piest <- mean(dat$Y[dat$X==1])*(1-plogis(deltagrid+betaest)) }
  if(x==1 & y==1){ piest <- mean(dat$Y[dat$X==1]==1) - mean(dat$Y[dat$X==0]==1)*plogis(deltagrid+betaest) }
  if(x==0 & y==0){ piest <- 1 - mean(dat$Y[dat$X==0]==1) - mean(dat$Y[dat$X==1]==0)*(1-plogis(deltagrid+betaest)) }
  
  U_q <- as.numeric(dat$X==1-x & dat$Y==y) - q
  U_r <- as.numeric(dat$X==1-x) - r
  U_S <- des*res
  
  U_beta <- as.numeric(dat$X==1-x & dat$Y==y)/q * predx - plogis(betaest)
  
  if(x==y){
    otherq <- mean(dat$X==x & dat$Y==y)
    U_otherq <- as.numeric(dat$X==x & dat$Y==y) - otherq
  }
  
  dU_q <- c(-1,rep(0,3+K))
  dU_r <- c(0,-1,rep(0,2+K))
  dU_S <- matrix(0,K,K+4)
  dU_S[,3:(3+K-1)] <- solve(-bread(fit))
  dU_beta <- c(mean(-1/q^2*as.numeric(dat$X==1-x & dat$Y==y)*predx),
               0,
               colMeans(as.numeric(dat$X==1-x & dat$Y==y)/q * model.matrix(glm(formula,data=newdatx,family="binomial")) * predx/(1+exp(predict(fit,newdatx))) ),
               -dplogis(betaest),
               0)
  
  if(x==1 & y==0){ 
    U_pi <- q/r * plogis(betaest + deltagrid) - piest
    dU_pi <- c(1/r * plogis(betaest+deltagrid),
               -q/r^2 * plogis(betaest+deltagrid),
               rep(0,K),
               q/r*plogis(betaest+deltagrid)/(1+exp(betaest+deltagrid)),
               -1)
  }
  if(x==0 & y==1){ 
    U_pi <- q/r * (1-plogis(betaest + deltagrid)) - piest
    dU_pi <- c(1/r * (1-plogis(betaest+deltagrid)),
               -q/r^2 * (1-plogis(betaest+deltagrid)),
               rep(0,K),
               -q/r*plogis(betaest+deltagrid)/(1+exp(betaest+deltagrid)),
               -1)
  }
  if(x==1 & y==1){ 
    U_pi <- otherq/(1-r) - q/r * plogis(betaest + deltagrid) - piest
    dU_pi <- c(-1/r * plogis(betaest+deltagrid),
               otherq/(1-r)^2 + q/r^2 * plogis(betaest+deltagrid),
               rep(0,K),
               -q/r*plogis(betaest+deltagrid)/(1+exp(betaest+deltagrid)),
               -1)
  }
  if(x==0 & y==0){ 
    U_pi <- otherq/(1-r) - q/r * (1-plogis(betaest + deltagrid)) - piest
    dU_pi <- c(-1/r * (1-plogis(betaest+deltagrid)),
               otherq/(1-r)^2 + q/r^2 * (1-plogis(betaest+deltagrid)),
               rep(0,K),
               q/r*plogis(betaest+deltagrid)/(1+exp(betaest+deltagrid)),
               -1)
  }
  
  if(x !=y ) { 
    U <- cbind(U_q,U_r,U_S,U_beta,U_pi)
    dU <- rbind(dU_q,dU_r,dU_S,dU_beta,dU_pi)
  }
  if(x ==y ) { 
    U <- cbind(U_q,U_r,U_otherq,U_S,U_beta,U_pi)
    dU_otherq <- rep(0,4+K)
    dU <- rbind(dU_q,dU_r,dU_otherq,dU_S,dU_beta,dU_pi)
    dU <- cbind(dU[,1:2],c(0,0,-1,rep(0,K),0,1/(1-r)),dU[,3:ncol(dU)]) 
  }
  
  cov <- solve(dU) %*% var(U) %*% t(solve(dU))/N
  
  list(betaest=betaest,
       piest=piest,
       U=U,
       dU=dU,
       cov=cov,
       varpi=cov[ncol(cov),ncol(cov)])
}
