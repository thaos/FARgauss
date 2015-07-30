packages <- c("boot", "quantreg", "evmix", "ismev", "parallel","foreach","doMC","snow")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	  install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.us.r-project.org")
}
library(boot)
library(quantreg)
library(evmix)
library(ismev)
library(parallel)
library(foreach)
library(doMC)
library(snow)

fit_gauss <- function(y,ydat,init=NULL){
	if (!is.null(init))
		mle=init
	else {
		tc=try({
			y.fit=lm(y~mua,data=ydat)
			fit_res=residuals(y.fit)
			var.fit=lm(fit_res^2~mua,data=ydat)
			print(y.fit)
		})
		if(class(tc)=="try-error" ){
			print("F********* n°1")
		}
		mle <- c(coefficients(y.fit),coefficients(var.fit))
	}
	gauss_lik <- function(mle){
		require(mvtnorm)
		mu0=mle[1]
		mu1=mle[2]
		sig0=mle[3]
		sig1=mle[4]
		mu=c(with(ydat,mu1*mua+mu0))
		sig2=with(ydat,diag(sig1*siga+sig0))
		-dmvnorm(y,mu,sig2,log=TRUE)
	}
	tc=try({
		y.fit2=nlminb(start=mle,gauss_lik)
	})
	if(class(tc)=="try-error" ){
		print("F********* n°2")
		browser()
	}
	y.fit2
}
	

getP <- function(p2,y.fit,ydat,to.plot=FALSE){
	xcord=p2[1]
	ycord=p2[2]
	x=ydat$year
	covariate=ydat$mua
	if (xcord < min(x) | xcord > max(x))
		stop(" Point outside of data range ")
	xcloser=which((abs(x-xcord))==min(abs(x-xcord)))
	xcloser=min(unique(x[xcloser]))
	ccloser=mean(unique(covariate[which(x==xcloser)]))
	mle <- y.fit$results$par
	mu=mle[1]+ccloser*mle[2]
	sigma2=mle[3]+ccloser*mle[4]
	shape=mle[5]
	1-evir::pgev(ycord,mu=mu,sigma=sigma,xi=shape)
	res=pnorm(ycord,mean=mu,sd=sqrt(sigma),lower.tail=FALSE)
	res=c(res,mu,sigma,shape)
	names(res)=c("p","mu","sigma","shape")
	res
}
# getP(c(1920,110),y.fit,ydat,to.plot=TRUE)

getFAR <- function(p1,x2,y.fit,ydat){
	nrob1=getP(p1,y.fit,ydat)
	p2=c(x2,p1[2])
	prob2=getP(p2,y.fit,ydat)
	if(prob1[1]==0 & prob2[1]==0)
		FAR=1
	else
	FAR=1-(prob2[1]/prob1[1])
	res=c(FAR,prob2,prob1)
	names(res)=c("FAR","p0","mu0","sc0","sh0","p1","mu1","sc1","sh1")
	res
}
# getFAR(c(1920,110),2102,y.fit,ydat)

getFAR.theo=function(xp,t0,t1,mu,sigma){
	p0=pnorm(xp,mean=mu[t0],sd=sigma[t0],lower.tail=FALSE)
	p1=pnorm(xp,mean=mu[t1],sd=sigma[t1],lower.tail=FALSE)
	print(p0/p1)
	if(p0==0 & p1==0)
		p0p1=1
	else
		p0p1=p0/p1
	res=1-p0p1
	attr(res,"p0")=p0
	attr(res,"p1")=p1
	res
}

FARBoot <- function(ydat,indice,p1,x2){
	y.fit.dat=fit_gauss(ydat$y,ydat)
	init=as.list(y.fit.dat$par)
	data.b=ydat[indice,]
	y.fit=fit_gauss(data.b$y,data.b,initial=init)
	getFAR(p1,x2,y.fit,ydat)
}


FARBoot.Spline <- function(ydat,indice,p1,x2){
	data.b=ydat[indice,]
	covariate=with(data.b,predict(rq( y~ bs(year, df=3,degree=3), tau=0.5)))
	data.b$mua=covariate
	data.b$siga=covariate
	y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga)
	getFAR(p1,x2,y.fit,ydat)
}

function2draw <- function(condition){
	function(ydat,nrepet,year){
		cond.e=eval(condition)
		sb=which(cond.e)
		sample(sb,size=nrepet,replace=TRUE)
	}
}


rbind_list <- function(liste){
	if(is.data.frame(liste)|length(liste)==1)
	   return(liste)
	if(length(liste)>1){
		liste[[2]]= rbind(liste[[1]],liste[[2]])
		return(rbind_list(liste[-1]))
	}
}

cond=quote(ydat$year==year)
draw_same_year <- function2draw(cond)
		
cond2=function(tol){
	substitute(ydat$year <= year+tol &ydat$year >= year-tol,list(tol=tol))
}
draw_around_year <- function2draw(cond2(5))

year2mu <- function(year,ydat){
	i=which(year==ydat$year)
	ydat[i[1],"mua"]
}
cond3 <- function(tol){
	substitute(ydat$mua <= year2mu(year,ydat)+tol &ydat$mua >= year2mu(year,ydat)-tol,list(tol=tol))
}
draw_around_mua <- function2draw(cond3(0.1))

FARBoot_gen  <- function(to_draw){
	function(ydat,indice,p1,x2){
		y.fit.dat=fevd(y,ydat,location.fun=~mua,scale.fun=~siga,method="MLE")
		init=as.list(y.fit.dat$results$par)
		years=aggregate(y~year,data=ydat,FUN=length)
		names(years)=c("year","eff")
		indice=c(with(years,mapply(to_draw,nrepet=eff,year=year,MoreArgs=list(ydat=ydat),SIMPLIFY=TRUE)))
		data.b=ydat[indice,]
		years.v=c(with(years,mapply(rep,x=year,times=eff)))
		data.b$year=years.v
		y.fit=fevd(y,data.b,location.fun=~mua,scale.fun=~siga,method="MLE",initial=init)
		getFAR(p1,x2,y.fit,ydat)
	}
}
FARBoot_ax=FARBoot_gen(draw_around_mua) 
FARBoot_sy=FARBoot_gen(draw_same_year) 
FARBoot_ay=FARBoot_gen(draw_around_year) 
