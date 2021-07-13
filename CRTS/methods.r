argv <- commandArgs(TRUE) # file, max.scale

# assume data is time-ordered
dat0 = as.matrix(read.table(argv[1],sep=" "))
#dat0 = as.matrix(read.table("1001004000475.dat",sep=" "))
obsx = dat0[,2]
n = length(obsx) - 2
times = dat0[,1]
times = times - times[1]
delta = diff(times)
mean.delta = (times[(n + 2)]) / (n + 1)

# results
results <- list()

### methods
### slepian wavelet filters (no approximation)
#library(waveslim, lib.loc="~mjg/lib/R/")
library(waveslim)
slepian.wf<-function(cons,jscale, mu, K=1, mesh)
{

## length of mesh should be equal to mm

  fl=1/2^(jscale+1)/mu
  fr=1/2^jscale/mu
  mm=(2^jscale)*cons

  if(mm!=length(mesh)) stop("check filter length")

  mat = matrix( rep(mesh,mm), mm,mm) - matrix( rep(mesh,mm), mm,mm, byrow=T)
  slep = (sin(2*pi*mat*fr)- sin(2*pi*mat*fl))/(pi*mat)
  diag(slep)=2*(fr-fl)
  slep= (diag(1,mm,mm)-1/mm)%*%slep%*%(diag(1,mm,mm)-1/mm)
  eigen(slep, symmetric=T)$vectors[seq(mm,1,-1),1:K]/sqrt(2^jscale*mu)
}

cons=1
max.scale = as.integer(argv[2])
sMj=sLj=rep(0, max.scale)
for(j in 1:max.scale) sLj[j] = 2^j*cons
for(j in 1:max.scale) sMj[j] = n-sLj[j] +1
mu=mean.delta

dpss <- vector("list", max.scale)
lamplus<-vector("list",max.scale)

for(jscale in 1:max.scale){  
  result <- tryCatch({
    dpss[[jscale]]<-dpss.taper(sMj[jscale],5,nw=3.5)
    lamplus[[jscale]]<-rep(1,sMj[jscale])%*%dpss[[jscale]]
    0
  }, warning = function(warn) {
     return(0)
   }, error = function(err) {
     return(-1)
   })
  if (result == -1) max.scale = max.scale - 1 # max.scale too big error
}

sw.var=sw.dis=tau=rep(0,max.scale)
wc=matrix(NA, max.scale,n)

for(jscale in 1:max.scale){

  result <- tryCatch({
    for(tpt in 1: sMj[jscale]){
      index=c((tpt+1):(tpt+sLj[jscale]))
      wf=slepian.wf(cons,jscale, mu, 1, times[index])	
      wc[jscale,tpt]=sum(obsx[index]*wf)
    }
    ser=wc[jscale, 1:sMj[jscale]]^2

    sw.var[jscale]=mean(ser)
   
    J= ser%*%dpss[[jscale]]
    u0= J%*%t(lamplus[[jscale]])/sum(lamplus[[jscale]]^2)
    sw.dis[jscale] = mean((J-c(u0)*lamplus[[jscale]])^2)/sMj[[jscale]]
    tau[jscale] = 2^(jscale - 1)*mean.delta
  }, warning = function(warn) {
    return(0)
  }, error = function(err) {
    return(-1)
  })
  if (result == -1) max.scale = max.scale - 1
}

#round( cbind(tau, sw.var, sw.dis), 9)
results['slep'] <- list(cbind(tau, sw.var, sw.dis))


### continuous time autoregressive model
suppressMessages(library(cts))
max.order = round(sqrt(n + 2))
fit <- vector("list", max.order)
aic = bic = 0
i = last.order = 1
while (i <= max.order) {
  scale = 2 * pi / mean.delta
  result = 1
  while (result > 0) {
    result <- tryCatch({
      fit[[i]] <- car(x = times, y = obsx, scale = scale, order = i, ctrl = car_control( vri = TRUE, ccv = "MNCT"))
      aic[i] <- fit[[i]]$aic
      bic[i] <- fit[[i]]$bic
      return(0)
    }, warning = function(warn) {
      return(0)
    }, error = function(err) {
      errcode = 0
      if (err$message == "ROOT WITH POSITIVE REAL PART") errcode = -1
      if (err$message == "program fails in rpoly") errcode = -1
      return(errcode)
    })
    if (result == -1) { # Increase scale by 10% to avoid positive real part root
      scale = scale * 1.1
      if (scale < 1) result = 1 # Exit condition
    }
    if (result < 0) {
      last.order = i
      i = max.order
    }
  }
  i = i + 1
}
order = which.min(bic)
results['car.basicfit'] <- list(fit[[1]]$b)
results['car.basicfit.scale'] <- 2 * pi / mean.delta
results['car.bestfit.orders'] <- list(c(order, last.order))
results['car.bestfit'] <- list(fit[[order]]$b)


### continuous time HMM with a single Gaussian per state
library(msm)
dofit <- function(data, order){
  qmat = matrix(runif(order * order), ncol = order)
  emFunc = hmmNorm(mean = 0.0, sd = 1.0)
  model = list(emFunc)[rep(1, order)]
  fit <- suppressWarnings(msm(V2 ~ V1, data = data, qmatrix = qmat, hmodel = model, hessian = FALSE))
}
data = as.data.frame(dat0)
data$V1 = data$V1 - data$V1[1]
data$V2 = data$V2 - mean(data$V2)
max.order = 10
res <- vector("list", max.order)
bic = 10000. # Some large value
for (i in 2:max.order) {
  res[[i]] <- try(dofit(data, i), TRUE)
#  res[[i]] <- suppressWarnings(dofit(data, i))
  bic[i] <- res[[i]]$minus2loglik + log(nrow(data)) * 2 * i
}
states = which.min(bic)
results['cthmm.order'] <- states
results['cthmm.fit'] <- list(res[[states]])
#print('CTHMM:')
#res[[states]])


### Teraesvarta nonlinearity
library(tseries)
fit <- terasvirta.test(times, obsx)
results['teraes.stat'] <- fit$statistic
results['teraes.pval'] <- fit$p.value


### Lyapunov exponent
freq = 1
N <- length(obsx)
Ly <- numeric(N - freq)
for (i in 1:(N - freq)) {
  idx <- order(abs(obsx[i] - obsx))
  idx <- idx[idx < (N-freq)]
  j <- idx[2] # Find nearest neighbour
  Ly[i] <- log(abs((obsx[i + freq] - obsx[j + freq]) / (obsx[i] - obsx[j]))) / freq
  if (is.na(Ly[i]) | Ly[i]==Inf | Ly[i]==-Inf)
    Ly[i] <- NA
}
Lyap <- mean(Ly,na.rm=TRUE)
fLyap <- exp(Lyap)/(1+exp(Lyap))
results['lyap'] <- fLyap


### Ashish's statistics
suppressMessages(library(e1071))

# convert matrix to dataframe
mat = data.frame(dat0)
names(mat)[1]<-paste("mjd")
names(mat)[2]<-paste("mag")
names(mat)[3]<-paste("magerr")

# Given a lightcurve, produce the statistics

lcstats <- function(olc){
  nobs <- nrow(olc)
  if(nobs < 5) return(NA)

  sorti <- order(olc$mjd)
  olc <- olc[sorti,]
  olc$measgrp <- factor(c(0,cumsum(ifelse(diff(olc$mjd) < 1, 0, 1))))

  # whole measures
  meds <- median(olc$mag)
  shov <- mean(abs(diff(olc$mag)))
  maxdiff <- max(abs(diff(olc$mag)))
  
  # fitted curve measures
  rmod <- lcfit(olc)
  dayrange <- max(olc$mjd) - min(olc$mjd)
  totvar <- sum(abs(diff(rmod$fitted)))/dayrange
  famp <- max(rmod$fitted) - min(rmod$fitted)
  
  # residual measures
  outl <- max(abs(rmod$residuals)/sd(rmod$residuals))
  skewres <- skewness(rmod$residuals)
  std <- sd(rmod$residuals)  
  
  # cluster measures
  if(length(levels(droplevels(olc$measgrp))) > 1){
    amod <- lm(mag ~ measgrp, olc)
    lsd <- log(summary(amod)$sigma)
    mdev <- max(abs(residuals(amod)))
  }else{
    lsd <- log(sd(olc$mag))
    mdev <- max(abs(olc$mag-mean(olc$mag)))
  }
  
  c(meds=meds, shov=shov, maxdiff=maxdiff, nobs=nobs, totvar=totvar, famp=famp, outl=outl, skewres=skewres, std=std, lsd=lsd, mdev=mdev)
}

# given a light curve, fit a curve and return fitted values and residuals

lcfit <- function(olc){
  if(nrow(olc) > 10){ # This is the Lowess version
    rmod <- tryCatch(loess(mag ~ mjd, olc),warning=function(x){return(lm(mag ~ mjd, olc))},error=function(x){return(lm(mag ~ mjd, olc))})
  }else{
    rmod <- lm(mag ~ mjd, olc)
  }
  rmod
}

results['aam'] <- list(lcstats(mat))


### Output results
print(results)
