# Define Log-Likelihood
log.like.fxn <- function(par,data,scaled.cov,y0){
  alpha <- par[["alpha"]] # Slope
  if(!is.null(y0)){y0<-y0}else{
    y0 <- par[["y0"]] # Intercept
  }
  rows <- nrow(data)
  log.like <- 0
  for(row in 1:rows){
    SBP_Change <- data$SBP_Change[[row]]/10
    COG_Change <- data$COG_Change[[row]]
    SBP_Change.SE <- data$SBP_Change.SE[[row]]/10
    COG_Change.SE <- data$COG_Change.SE[[row]]
    #ss <- data$Sample.Size[[row]]
    #if(SBP_Change.SE==0){log.like <- log.like}else{
    years <- data$Intervention.Length.Months[[row]]/12
    predicted.cog <- (alpha*SBP_Change+y0)*years
    se.predicted.cog <- abs(SBP_Change.SE*alpha*years)
    cov <- scaled.cov*COG_Change.SE*se.predicted.cog
    err <- sqrt(COG_Change.SE^2+se.predicted.cog^2-
                  2*cov)
    deg.free <- Inf
    log.prob <- dt((COG_Change-predicted.cog)/
                     err,
                   df=deg.free,log=T)
    log.like <- log.like + log.prob 
    #}
  }
  if(is.finite(log.like))log.like
  else{-10000}
}