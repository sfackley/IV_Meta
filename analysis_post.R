rm(list=ls())
require('knitr')
require('tidyr')
require('dplyr')
require('gtools')
require('plotly')
require('gdata')
require('ggplot2')
require('magrittr')
require('reshape')
require('rlang')
require('pracma')
require('ggthemes')
require('binom')
require('scales')
require('cowplot')
require('grid')
load("data.RData")

dat <- data.frame(treatment=c("Donanemab","Placebo"))
dat$suvr.24 <- c(-67.83-2,-2)/100
dat$suvr.52 <- c(-82.30,0)/100
dat$suvr.76 <- c(-85.06-1,-1)/100
dat$se.suvr.24 <- c(2.234457,2.234457)/100
dat$se.suvr.52 <- c(2.411234,2.411234)/100
dat$se.suvr.76 <- c(2.736503,2.736503)/100
dat$cog.24 <- c(-0.71,-0.7)
dat$cog.52 <- c(-1.17,-1.56)
dat$cog.76 <- c(-2.35,-2.98)
dat$se.cog.24 <- c(0.292,0.295)
dat$se.cog.52 <- c(0.321,0.321)
dat$se.cog.76 <- c(0.386,0.390)
dat$suvr.total.24 <- c(115,111) 
dat$suvr.total.52 <- c(92,91)
dat$suvr.total.76 <- c(90,91)
dat$cog.total.24 <- c(115,111)
dat$cog.total.52 <- c(92,91) 
dat$cog.total.76 <- c(90,91) 

## Reformat 
select <- dplyr::select
extract <- tidyr::extract
rename <- dplyr::rename
donanemab.dat <- dat %>% 
  select("treatment",everything())
donanemab.dat %<>%
  gather(value,key,suvr.24:cog.total.76) %>% 
  extract("value",
          c("name","week"),
          "(.+).([[:digit:]]{2})") %>%
  spread(name,key) 
save(donanemab.dat,file="data_donanemab.Rdata")


## DEFINE LIKELIHOOD FUNCTION
## See https://en.wikipedia.org/wiki/Welch's_t-test
## In statistics, Welch's t-test, or unequal variances t-test, 
## is a two-sample location test which is used to test the
## hypothesis that two populations have equal means.
log.like.fxn <- function(par,data){
  alpha <- par[["alpha"]]
  y0 <- par[["y0"]]
  rows <- nrow(data)
  log.like <- 0
  for(row in 1:rows){
    suvr <- data$suvr[[row]]
    cog <- data$cog[[row]]
    se.cog <- data$se.cog[[row]]
    se.suvr <- data$se.suvr[[row]]
    cog.total <- data$cog.total[[row]]
    suvr.total <- data$suvr.total[[row]]
    predicted.cog <- alpha*suvr+y0
    se.predicted.cog <- se.suvr*alpha
    err <- sqrt(se.cog^2+se.predicted.cog^2)
    deg.free <- err^4/(se.cog^4/(cog.total-1)+
                         se.predicted.cog^4/(suvr.total-1))
    #print(deg.free)
    y.prob <- dt((cog-predicted.cog)/
                   err,
                 df=deg.free)
    row.prob <- y.prob
    log.like <- log.like + log(row.prob)
  }
  if(is.finite(log.like))log.like
  else{-10000}
}

## FUNCTION TO GET PARAMETERS
get.pars <- function(data,adas.data=F){
  par0 <- lm(cog~suvr,data=data)$coefficients
  par0 <- rev(par0)
  names(par0) <- c("alpha","y0")
  log.like.fxn(par0,data)
  optim.fun <- function(x)(-log.like.fxn(x,data=data))
  if(adas.data){
    lower.y0 <- -20
    upper.y0 <- 20
    lower.alpha <- -500
    upper.alpha <- 200
  }else{
    lower.y0 <- -10
    upper.y0 <- 10
    lower.alpha <- -60
    upper.alpha <- 60    
  }
  optim.out <- optim(par0,optim.fun,
                     lower=c(lower.alpha,lower.y0),
                     upper=c(upper.alpha,upper.y0),
                     method="L-BFGS-B")
  optim.out$convergence
  pars <- optim.out$par
  hess <- hessian(optim.fun,pars)
  fisher.info <- solve(hess)
  sigma <- sqrt(diag(fisher.info))
  z.crit <- qnorm(0.975)
  lower <- pars - z.crit*sigma
  upper <- pars + z.crit*sigma
  pars.table <- list(pars=pars,
                     lower=lower,
                     upper=upper)
  pars.table
}

# One trial of bapineuzumab only gives diff in cog
# relative to placebo.
# So, I constrain the intercept to be zero.
# bap.data1
log.like.fxn.bap <- function(alpha,data){
  # -intercept = alpha*bap.data1$suvr[2]
  pars <- c("alpha"=alpha,"y0"=-alpha*bap.data1$suvr[2])
  log.like.fxn(pars,data)
}

get.pars.bap <- function(data=bap.data1){
  par0 <- lm(cog~suvr,data=data)$coefficients
  par0 <- rev(par0)
  names(par0) <- c("alpha","y0")
  par0 <- par0[[1]]
  log.like.fxn.bap(par0,data)
  optim.fun <- function(x)(-log.like.fxn.bap(x,data=data))
  optim.out <- optim(par0,optim.fun,
                     lower=c(-60),
                     upper=c(50),
                     method="L-BFGS-B")
  optim.out$convergence
  pars <- optim.out$par
  hess <- hessian(optim.fun,pars)
  fisher.info <- solve(hess)
  sigma <- sqrt(diag(fisher.info))
  z.crit <- qnorm(0.975)
  lower <- pars - z.crit*sigma
  upper <- pars + z.crit*sigma
  pars.table <- list(pars=pars,
                     lower=lower,
                     upper=upper)
  pars.table
}

## BAN2401

log.like.fxn.ban <- function(par,data=ban.dat){
  alpha <- par[["alpha"]]
  beta <- par[["beta"]]
  y0 <- par[["y0"]]
  rows <- nrow(data)
  log.like <- 0
  for(row in 1:rows){
    suvr <- data$suvr[[row]]
    cog <- data$cog[[row]]
    se.cog <- data$se.cog[[row]]
    se.suvr <- data$se.suvr[[row]]
    cog.total <- data$cog.total[[row]]
    suvr.total <- data$suvr.total[[row]]
    prop.e4 <- data$prop.e4.pos[[row]]
    month <- data$month[[row]]
    month <- as.numeric(month)
    # Adding the effect of time for APOE-E4
    # 12 months is actually 53 weeks
    # 18 months is actually 79 weeks
    predicted.cog <- alpha*suvr+y0+
      beta*prop.e4*
      (53/52)^(month==12)*
      (79/52)^(month==18)*
      (78/52)^(month==1)
    se.predicted.cog <- se.suvr*alpha
    err <- sqrt(se.cog^2+se.predicted.cog^2)
    deg.free <- err^4/(se.cog^4/(cog.total-1)+se.predicted.cog^4/(suvr.total-1))
    y.prob <- dt((cog-predicted.cog)/
                   err,
                 df=deg.free)
    row.prob <- y.prob
    log.like <- log.like + log(row.prob)
  }
  if(is.finite(log.like))log.like
  else{-100000000}
}

get.pars.ban <- function(data=ban.dat){
  par0 <- lm(cog~suvr+prop.e4.pos,data=data)$coefficients
  par0 <- rev(par0)
  names(par0) <- c("beta","alpha","y0")
  log.like.fxn.ban(par0,data)
  optim.fun <- function(x)(-log.like.fxn.ban(x,data=data))
  optim.out <- optim(par0,optim.fun,
                     lower=c(-10,-60,-10),
                     upper=c(10,50,10),
                     method="L-BFGS-B")
  optim.out$convergence
  pars <- optim.out$par
  print(optim.out)
  hess <- hessian(optim.fun,pars)
  print(hess)
  fisher.info <- solve(hess)
  sigma <- sqrt(diag(fisher.info))
  z.crit <- qnorm(0.975)
  lower <- pars - z.crit*sigma
  upper <- pars + z.crit*sigma
  pars.table <- list(pars=pars,
                     lower=lower,
                     upper=upper)
  pars.table
}

## Get parameter estimates
bex.pars <- get.pars(bex.data)
sol.pars <- get.pars(sol.data)
sol.pars2 <- get.pars(sol.data2)
sol.pars3 <- get.pars(sol.data3)
ly.pars <- get.pars(ly.data)
ly.pars.mmse <- get.pars(ly.data.mmse)
ly.pars2 <- get.pars(ly2.data)
gan.pars <- get.pars(gan.data)
bap1.pars <- get.pars.bap(bap.data1) 
bap.e4.pars <- get.pars(bap.reg.data %>% filter(E4=="e4"))
bap.noe4.pars <- get.pars(bap.reg.data %>% filter(E4=="noe4"))
ver.pars <- get.pars(ver.data)
ver.pars2 <- get.pars(ver.data2)
ver.pars2.mmse <- get.pars(ver.data2.mmse)
ban.pars <- get.pars.ban(ban.dat)
ban.pars.mmse <- get.pars.ban(ban.dat.mmse)
ban.est <- ban.pars %>% unlist %>% t %>% 
  as.data.frame %>% select(contains("alpha"))
emerge.pars <- get.pars(emerge.data)
engage.pars <- get.pars(engage.data)
donanemab.pars <- get.pars(donanemab.dat)

log.like.all.data <- function(pars,
                                     bapdata1=bap.data1,
                                     bapdata2=bap.reg.data %>% filter(E4=="e4"),
                                     bapdata3=bap.reg.data %>% filter(E4=="noe4"),
                                     soldata1=sol.data,
                                     soldata2=sol.data2,
                                     soldata3=sol.data3,
                                     bexdata=bex.data,
                                     gandata=gan.data,
                                     lydata=ly.data.mmse,
                                     ly2data=ly2.data, 
                                     verdata=ver.data,
                                     verdata2=ver.data2.mmse,
                                     bandata=ban.dat.mmse,
                                     emergedata=emerge.data,
                                     engagedata=engage.data
){
  alpha=pars[["alpha"]]
  log.like <- log.like.fxn.bap(alpha,bapdata1) 
  #print(log.like)
  tmp.pars <- c(alpha=pars[["alpha"]],y0=pars[["y0bap2"]])
  log.like <- log.like + log.like.fxn(tmp.pars,bapdata2) 
  #print(log.like)
  tmp.pars <- c(alpha=pars[["alpha"]],y0=pars[["y0bap3"]])
  log.like <- log.like + log.like.fxn(tmp.pars,bapdata3)
  #print(log.like)
  #tmp.pars <- c(y0=pars[["y0sol1"]],alpha=pars[["alpha"]])
  #REDUNDANT log.like <- log.like + log.like.fxn(tmp.pars,soldata1) 
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0sol2"]],alpha=pars[["alpha"]])
  log.like <- log.like + log.like.fxn(tmp.pars,soldata2) 
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0sol3"]],alpha=pars[["alpha"]])
  log.like <- log.like + log.like.fxn(tmp.pars,soldata3)
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0bex"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,bexdata)
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0gan"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,gandata)
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0ly2"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,ly2data)
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0ver"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,verdata)
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0ver2"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,verdata2)
  #print(log.like)
  tmp.pars <- c(y0=pars[["y0ly"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,lydata)
  #
  tmp.pars <- c(alpha=pars[["alpha"]],beta=pars[["beta"]],y0=pars[["y0ban"]])
  log.like <- log.like + log.like.fxn.ban(tmp.pars,bandata) 
  #
  tmp.pars <- c(y0=pars[["y0em"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,emergedata)
  #
  tmp.pars <- c(y0=pars[["y0en"]],alpha=pars[["alpha"]]) 
  log.like <- log.like + log.like.fxn(tmp.pars,engagedata)
  #
  log.like
}
par0 <- c(alpha=-1.8886057,
          y0bap2=-4.5049627, y0bap3=-3.8224559,
          # y0sol1=-2.2844224,
          y0sol2=-7.8469647,
          y0sol3=-3.4214022,
          y0bex=1.0443271,y0gan=-2.9566891,y0ly2=-3.4493283,
          y0ver=-3.9046654,y0ver2=-3,
          y0ly=-1.5068869,y0ban=-2.0520561,
          y0em=-3.382042,y0en=-3.3971117,
          beta=-0.8594328)
log.like.all.data(par0)

optim.fun <- function(x)-log.like.all.data(x)
optim.out <- optim(par0,optim.fun,
                   lower=par0-10,
                   upper=par0+10,
                   method="L-BFGS-B",
                   control=list(factr=1e-8))
optim.out$convergence
pars <- optim.out$par
pars
hess <- hessian(optim.fun,pars)
fisher.info <- solve(hess)
sigma <- sqrt(diag(fisher.info))
z.crit <- qnorm(0.975)
lower <- pars - z.crit*sigma
upper <- pars + z.crit*sigma
all.pars <- list(pars=pars,lower=lower,upper=upper)

