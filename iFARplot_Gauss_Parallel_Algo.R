packages <- c("boot", "quantreg", "evmix", "ismev", "gWidgetstcltk", "gWidgets","evir","extRemes")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
	  install.packages(setdiff(packages, rownames(installed.packages())),repos="http://cran.us.r-project.org")
}
library(boot)
library(quantreg)
library(evmix)
library(ismev)
library(evir)
library(extRemes)
require(gWidgetstcltk) #use the gWidgetstcltk package from CRAN
require(gWidgets) #use the gWidgetstcltk package from CRAN
source("FAR_GAUSS_Algo.R")
source("Gauss_ICratio_mu_alabama_IR2.R")


iFARplot <- function() {
	thisEnv <- environment()
	FarHat=NULL
	repet=4
	years=1850:2200
	n=length(unique(years))
	changement=100
	t=seq(1,n)
	mu=(t>changement)*.020*(t-changement)+100
	sigma2=(t>changement)*.001*(t-changement)+1
	shape=rep(-0.10,n)
	mu=rep(mu,repet)
	sigma=sqrt(rep(sigma,repet))
	t=rep(t,repet)
	shape=rep(shape,repet)
	year=rep(years,repet)
	covariate=(mapply(qnorm,mean=mu,sd=sigma,MoreArgs=list("p"=0.5)))
	#On ajoute un bruit insignifiant pour éviter que la covariable soit exactementy identique sur la première période, sans quoi, ce la fait bugguer le profile de vraisemblance: division par zéro pour une différence en covariable nulle.
	covariate=covariate+rnorm(n*repet,sd=0.000001)
	covariate=scale(covariate,scale=FALSE)
	y=mapply(rnorm,1,mu,sigma)
	ydat=data.frame(year,y,mu,sigma,shape,rep(1,n*repet),covariate,rep(1,n*repet),covariate,rep(1,n*repet))
	names(ydat)=c("year","y","mu","sigma","shape","mub","mua","sigb","siga","xi")
	df.plot=ydat[order(t),]
	getDensity.theo=function(tx,mu,sigma,xi){
		ymin=min(y)
		ymax=max(y)
		x=seq(ymin,ymax,length.out=100)
		res=sapply(x,devd,loc=mu[tx],scale=sigma[tx],shape=xi[tx])
		res=cbind(x,res)
		res
	}
	getDensity=function(mu,sigma,xi){
		ymin=min(y)
		ymax=max(y)
		x=seq(ymin,ymax,length.out=100)
		res=sapply(x,devd,loc=mu,scale=sigma,shape=xi)
		res=cbind(x,res)
		res
	}
	updatePlot <- function(h,...){
		layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
		with(df.plot,plot(year,y,main="Synthetic Time Serie of Maxima"))
		startAdjust[]<<-seq(min(year),svalue(endAdjust),by=1)
		endAdjust[]<<-seq(svalue(startAdjust),max(year),by=1)
		minEvent=with(df.plot,min(y[year>svalue(startAdjust)&year<svalue(endAdjust)]))
		maxEvent=with(df.plot,max(y[year>svalue(startAdjust)&year<svalue(endAdjust)]))
		EventAdjust[] <<- seq(minEvent,maxEvent,length=100)
		abline(v=svalue(startAdjust),col="blue")
		text(x=svalue(startAdjust),y=109,pos=2,labels="t0",col="red")
		abline(v=svalue(endAdjust),col="blue")
		text(x=svalue(endAdjust),y=109,pos=4,labels="t1",col="red")
		abline(h=svalue(EventAdjust),col="blue")
		text(y=svalue(EventAdjust),x=2200,pos=1,labels="xp",col="red")
		xp=svalue(EventAdjust)
		t0=svalue(startAdjust)
		t1=svalue(endAdjust)
		i0=min(which((abs(df.plot$year-t0))==min(abs(df.plot$year-t0))))
		i1=min(which((abs(df.plot$year-t1))==min(abs(df.plot$year-t1))))
		d0=getDensity.theo(tx=i0,mu=df.plot$mu,sigma=df.plot$sigma,xi=df.plot$shape)
		d1=getDensity.theo(tx=i1,mu=df.plot$mu,sigma=df.plot$sigma,xi=df.plot$shape)
		cat("JeanClaudeRichard")
		if(!is.null(get("FarHat",envir=thisEnv))){
			cat("PierreMarieSolange")
			d0.prof=with(FarHat,getDensity(mu=Prof[2],sigma=Prof[3],xi=Prof[4]))
			d1.prof=with(FarHat,getDensity(mu=Prof[6],sigma=Prof[7],xi=Prof[8]))
			d0.boot=with(FarHat,getDensity(mu=Boot[2],sigma=Boot[3],xi=Boot[4]))
			d1.boot=with(FarHat,getDensity(mu=Boot[6],sigma=Boot[7],xi=Boot[8]))
			print(all(d0.boot==d0.prof))
		}	
		r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=df.plot$mu,sigma=df.plot$sigma,xi=df.plot$shape)
		p0.theo=attr(r.theo,"p0")
		p1.theo=attr(r.theo,"p1")
		plot(d0[,1],d0[,2],main="Theoric density at time t0",type="l",xlab="",ylab="density",ylim=c(0,0.5),sub=paste("theoric p0:" ,format(p0.theo,digits=2),sep=""))
		if(!is.null(get("FarHat",envir=thisEnv))){
			lines(d0.prof[,1],d0.prof[,2],col="blue")
			lines(d0.boot[,1],d0.boot[,2],col="green")
			legend("topleft",legend=c("Profile","Boot"),lty=1,lwd=2,col=c("blue","green"),bty="n")
		}
		abline(v=xp)
		text(x=xp,y=0.45,pos=4,labels="xp",col="red")
		plot(d1[,1],d1[,2],main="Theoric density at time t1",type="l",xlab="",ylab="density",ylim=c(0,0.5),sub=paste("theoric p1:" ,format(p1.theo,digits=2),sep=""))
		if(!is.null(get("FarHat",envir=thisEnv))){
			lines(d1.prof[,1],d1.prof[,2],col="blue")
			lines(d1.boot[,1],d1.boot[,2],col="green")
			legend("topleft",legend=c("Profile","Boot"),lty=1,lwd=2,col=c("blue","green"),bty="n")
		}
		abline(v=xp)
		text(x=xp,y=0.45,pos=4,labels="xp",col="red")
	}
	ComputeFAR <- function(h,...){
		svalue(resText) <<- paste("Computing...\n",sep="")
		xp=svalue(EventAdjust)
		t0=svalue(startAdjust)
		t1=svalue(endAdjust)
		i0=min(which((abs(df.plot$year-t0))==min(abs(df.plot$year-t0))))
		i1=min(which((abs(df.plot$year-t1))==min(abs(df.plot$year-t1))))
		print(paste(t0,i0))
		print(paste(t1,i1))
		y.fit=fevd(y,df.plot,location.fun=~mua,scale.fun=~siga)
		cat("y fitted")
		r.ic=gauss.ratio.ic.mu(xp=xp,t0=i0,t1=i1,ydat=df.plot)
		r.param=r.ic[4:length(r.ic)]
		r.ic[1:3]=1-r.ic[1:3]
		r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=df.plot$mu,sigma=df.plot$sigma)
		p0.theo=attr(r.theo,"p0")
		p1.theo=attr(r.theo,"p1")
		svalue(resText) <<- paste("FAR pour t0=",format(t0,digits=2)," et t1=",format(t1,digits=2)," et seuil=",format(xp,digits=2),"\n\n",sep="")
		insert(resText,paste("theoric FAR:\n", format(r.theo,digits=2),"\n",sep=""))
		insert(resText,paste("theoric p0\t\ttheoric p1\n", format(p0.theo,digits=2),"\t\t",format(p1.theo,digits=2),"\n",sep=""))
		prof.txt=paste("Profile Likelihood FAR:\nEstimate\t\tIC95_inf\t\tIC95_sup\n",format(r.ic[2],digits=2),"\t\t",format(r.ic[3],digits=2),"\t\t",format(r.ic[1],digits=2),"\n",sep="")
		insert(resText,prof.txt)
		insert(resText,paste("p0.esti\t\tp1.esti\n", format(r.ic[4],digits=2),"\t\t",format(r.ic[8],digits=2),"\n",sep=""))
		cat("Computing Bootstrap IC\n")
		boot.res=boot(data=df.plot,statistic=FARBoot,R=250,p1=c(t1,xp),x2=t0,parallel="snow")
		r.boot=colMeans(boot.res$t)
		if(dev.cur()== dev.next())
			dev.new()
		else
			dev.set(dev.next())
		par(mfrow=c(1,3))
		hist(boot.res$t[,1],main="BootSample : FAR",xlab="FAR",breaks="FD",freq=FALSE)
		abline(v=r.theo,col="red")
		abline(v=r.boot[1],col="green")
		hist(boot.res$t[,2],main="BootSample : p0",xlab="p0",breaks="FD",freq=FALSE)
		abline(v=p0.theo,col="red")
		abline(v=r.boot[2],col="green")
		hist(boot.res$t[,6],main="BootSample : p1",xlab="p1")
		#                 hist(boot.res$t[,6],main="BootSample : p1",xlab="p1",breaks="FD",freq=FALSE)
		abline(v=p1.theo,col="red")
		abline(v=r.boot[6],col="green")
		dev.set(dev.prev()) 
		alpha=0.05
		ic.boot=apply(boot.res$t,2,quantile,p=c(alpha/2,1-alpha/2))
		i0=min(which((abs(df.plot$year-t0))==min(abs(df.plot$year-t0))))
		i1=min(which((abs(df.plot$year-t1))==min(abs(df.plot$year-t1))))
		r.theo=getFAR.theo(xp=xp,t0=i0,t1=i1,mu=df.plot$mu,sigma=df.plot$sigma)
		r.theo=c(r.theo,attr(r.theo,"p0"),attr(r.theo,"p1"))
		boot.ic=as.data.frame(t(rbind(r.boot,ic.boot)[,c(1,2,6)]))
		boot.ic=cbind(r.theo,boot.ic)
		names(boot.ic)=c("Theoric","Estimate","LowerIC","UpperIC")
		rownames(boot.ic)=c("FAR","p0","p1")
		print(boot.ic)
		boot.txt=paste("Resampling Bootstrap FAR:\nEstimate\t\tIC95_inf\t\tIC95_sup\n",format(boot.ic[1,2],digits=2),"\t\t",format(boot.ic[1,3],digits=2),"\t\t",format(boot.ic[1,4],digits=2),"\n",sep="")
		insert(resText,boot.txt)
		insert(resText,paste("p0.esti\t\tp1.esti\n",format(boot.ic[2,2],digits=2),"\t\t",format(boot.ic[3,2],digits=2),"\n",sep=""))
		assign("FarHat",list("Prof"=r.param,"Boot"=r.boot[2:length(r.boot)]),envir=thisEnv)
		updatePlot()
		assign("FarHat",NULL,envir=thisEnv)
	}
	##The widgets
	win <- gwindow("Fraction of Attributable Risk")
	gp <- ggroup(horizontal=FALSE, cont=win)
	tmp <- gframe("Results", container=win, expand=TRUE)
	resText <- gtext("FAR",cont=tmp,width=320,height=340)
	tmp <- gframe("Start-End", container=gp, expand=TRUE)
	startAdjust <- gslider(from=min(year),to=max(year),by=1, value=min(year),
				   cont=tmp, expand=TRUE,
				   handler=updatePlot)
	endAdjust <- gslider(from=min(year),to=max(year),by=1, value=max(year),
				   cont=tmp, expand=TRUE,
				   handler=updatePlot)
	minEvent=min(df.plot$y)
	maxEvent=max(df.plot$y)
	tmp <- gframe("event.threshold", container=gp, expand=TRUE)
	EventAdjust <- gslider(from=minEvent,to=maxEvent,length=100, value=minEvent,
				   cont=tmp, expand=TRUE,
				   handler=updatePlot)
	gbutton("ComputeFAR",container=tmp,handler=ComputeFAR)
	updatePlot()
}
