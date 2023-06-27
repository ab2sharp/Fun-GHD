library(Bessel)
library(stats)
library(MASS)
library(mvtnorm)
library(ghyp)
library(numDeriv)
library(mixture)








ddghyp <- function(x=NULL, par=NULL, logd=FALSE) {
# x  is a n * p matrix
	mu    = par$mu
	beta = par$beta

	if (length(mu) != length(beta) ) stop("mu and beta do not have the same length")
  omega1 = par$ol1[1]; lambda1 = par$ol1[2]; 
  omega2 = par$ol2[1]; lambda2 = par$ol2[2]; 

	gam = par$gam; psi = par$psi; eta = par$eta; d = length(psi); p = nrow(gam); n = nrow(x);
	
	gam.x = x %*% gam
	gam.xmu = sweep(gam.x,2,mu)
	
	gam.xmu2.ipsi = sweep( gam.xmu^2, 2, 1/psi, "*") 
	quad1 = apply(gam.xmu2.ipsi,1,sum)
 	quad2 = (apply(x^2, 1, sum)  - apply(gam.x^2,1,sum))/eta
 	
	pa1 = omega1 + sum( beta^2/psi )
	mx1 = omega1 + quad1

	pa2 = omega2 + 0
	mx2 = omega2 + quad2
	
	kx1 = sqrt(mx1*pa1); kx2 = sqrt(mx2*pa2)

	lvx1 = matrix(0, nrow=nrow(x), 3)
	lvx1[,1] = (lambda1 - d/2)*log(kx1)
	lvx1[,2] = log(besselK( kx1, nu=lambda1-d/2, expon.scaled = TRUE)) -kx1
	lvx1[,3] = as.numeric( gam.xmu %*% diag(1/psi) %*% beta )

	lv1 = numeric(3)
	lv1[1] = -1/2*sum(log(psi)) - d/2*(log(2)+log(pi)) 
	lv1[2] = omega1 - log(besselK( omega1, nu=lambda1, expon.scaled = TRUE)) 
	lv1[3] = (d/2 - lambda1)*log( pa1 )

	lvx2 = matrix(0, nrow=nrow(x), 3)
	lvx2[,1] = (lambda2 - (p-d)/2)*log(kx2)
	lvx2[,2] = log(besselK( kx2, nu=lambda2-(p-d)/2, expon.scaled = TRUE)) -kx2
	lvx2[,3] = 0 
	
	lv2 = numeric(3)
	lv2[1] = -1/2*(p-d)*log(eta) - (p-d)/2*(log(2)+log(pi))
	lv2[2] = omega2 - log(besselK( omega2, nu=lambda2, expon.scaled = TRUE)) 
	lv2[3] = ( (p-d)/2 - lambda2)*log( pa2 )
	
	val = apply(lvx1,1,sum) + apply(lvx2,1,sum)  + sum(lv1) + sum(lv2)
	if (!logd) val = exp( val )

	return(val)
}




llik <- function(data,gpar, delta=0) {
	logz = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
	for (k in 1:length(gpar$pi) ) logz[,k] = ddghyp(x=data, par=gpar[[k]], logd=TRUE)
	val = sum(log(apply(logz,1,function(z,wt=NULL) {
		return(sum(exp(z)*wt))
	},wt=gpar$pi)))

	return(val)
	}

MAP <- function(data, gpar, label=NULL) {
	w = weights(data=data, gpar=gpar, v=1)
	if (!is.null(label)) w = combinewk(weights=w, label= label)
	z = apply(w, 1, function(z) { z=(1:length(z))[z==max(z)]; return(z[1]) })
	z = as.numeric(z)
	return( z)	
	}


weights <- function(data=NULL, gpar=NULL, v=1) {
  G = length(gpar$pi)	
  if (G > 1) {
    zlog = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:G ) zlog[,k] =  ddghyp(x=data, par=gpar[[k]], logd=TRUE)
    w = t(apply(zlog, 1, function(z,wt,v) { 
      x=exp(v*(z + log(wt) -max(z)) );
      
     # x[x< 0.0001] = 0
      x =  x/sum(x)
      return( x ) 
    }, wt=gpar$pi,v=v ))
  } else w = matrix(1,nrow=nrow(data), ncol=G)
  return(w)
}
	
combinewk <- function(weights=NULL, label=NULL)	{
	# known is a numeric with 
	# 0 if unknown group membership 
	# 1,2,3,.. for label of known group
	if (is.null(label)) stop('label is null')
	kw     = label !=0
 	for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
	return(weights)	
}






ipar <- function(data, wt=NULL, d=NULL) {
	if (is.null(wt)) wt = rep(1,nrow(data))
	#Dimension of the data (dimensions are columns, observations are rows)
	p = ncol(data)
	val = list()
	#Generate random mean and skewness vectors
	val$mu    = rnorm(d, apply(data,2,weighted.mean, w=wt), sqrt(1/nrow(data)) )
	val$beta = rnorm(d, 0,  1/nrow(data) )

	temp = eigen(cov.wt(data, wt = wt, method="ML")$cov )
	if (d == 0) val$gam = NULL
	else if (d == 1) val$gam = matrix( temp$vectors[,1:d], p, 1)
	else val$gam = temp$vectors[,1:d]

	val$psi  = temp$values[1:d]
	val$eta = mean(temp$values[-c(1:d)])*10
		
	val$ol1  = c( 0.1, -p)
	val$ol2 = c( 0.1, -p)

	return(val)
}






#Randomly generated params -- calls ipar
rgpar <- function(data, g=2, w=NULL, d=NULL ) {
	if (is.null(w) ) w = matrix(1/g, nrow=nrow(data), ncol=g)
	val = list()
	for (k in 1:g) val[[k]] = ipar(data=data, wt=w[,k], d=d[k])
	val$pi = rep(1/g,g)
	val$d  = d
	return(val)
}	
	




igpar <- function(
	data=NULL, g=NULL, nstart=5, iter=10, mtol= NULL, mmax= NULL, 
	covtype= "VVVV", skewness=TRUE, dg=NULL 
) {
	z = numeric(nstart+3)	
	gparn <- list()
	
	gparn[[1]] = gpcmgpar(data=data, G=g, n= iter,  mtol=mtol, mmax= mmax, covtype= covtype,v=0, skewness= skewness, dg=dg )
	z[1] = llik(data, gparn[[1]])
	gparn[[2]] = kmeansgpar(data=data, G=g, n= iter, mtol=mtol, mmax= mmax, covtype= covtype, v=0, skewness= skewness, dg=dg )
	z[2] = llik(data, gparn[[2]])
	gparn[[3]] = hclustgpar(data=data, G=g, n= iter, mtol=mtol, mmax= mmax, covtype= covtype, v=0, skewness= skewness, dg=dg )
	z[3] = llik(data, gparn[[3]])

	for (i in 1:nstart) {
		#print(i)
		gpari = rgpar(data=data, g=g, d=dg )
		gpari = EMn(data=data, gpar0=gpari, G=g, mtol=mtol, mmax= mmax, covtype= covtype, skewness= skewness)$gpar
		gparn[[i+3]] = gpari
		z[i+3] = llik(data, gparn[[i+3]])
	}

	remove = is.infinite(z)
	z = z[!remove]
	gparn = gparn[!remove]
	gpar = gparn[[ order(z,decreasing=TRUE)[1] ]]	

	return(gpar)
}

	

gpcmgpar <- function(
	data=NULL, G=2, v=1, label=NULL, n=50, mtol= NULL, mmax= NULL, 
	covtype= NULL, skewness= FALSE, dg=NULL 
) {
	w = gpcm(data, G=G, mnames=NULL)$z
	gpar  = rgpar(data=data, g=G, w=w, d=dg)
	
	
	# do n iterations under these weights
	for (j in 1:n) try({ gpar = EMgrstep(data=data, gpar=gpar, v=v, label=label, w=w,  mtol=mtol, mmax= mmax, covtype= covtype, skewness= skewness) }, TRUE)
	try({ gpar = EMn(data=data, gpar0=gpar, label=label, n=n,  mtol=mtol,mmax= mmax, covtype= covtype, skewness= skewness )$gpar })
	return(gpar)
}

	
	
hclustgpar <- function(
	data=NULL, G=NULL, n=50, label=NULL, v=1,  mtol=NULL, mmax= NULL, covtype= NULL, skewness= FALSE, dg=NULL  
) {
	lw = cutree(hclust(dist(data), "complete"),k=G)
	z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
	gpar  = rgpar(data=data, g=G, w=z, d=dg)

	for (j in 1:n) try({ gpar = EMgrstep(data=data, gpar=gpar, v=v, label= label, w=z,  mtol=mtol, mmax= mmax, covtype= covtype, skewness= skewness )}, TRUE)
	try({ gpar = EMn(data=data, gpar=gpar, label= label, mtol=mtol, mmax= mmax, covtype= covtype, skewness= skewness )$gpar }, TRUE)

	return(gpar)
}
	
	
kmeansgpar <- function(
	data=NULL, G=NULL, n=10, label=NULL, v=1,  mtol= NULL, mmax= NULL, covtype= NULL, skewness= FALSE, dg=NULL  
) {
	#print("igpar3 1")
	lw = kmeans(data, centers=G, iter.max=100, nstart=10)$cluster
	z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
	gpar  = rgpar(data=data, g=G, w=z, d=dg)
	
	for (j in 1:n) try({ gpar = EMgrstep(data=data, gpar=gpar, v=v, label= label, w=z, mtol=mtol, mmax= mmax, covtype= covtype, skewness= skewness )}, TRUE)
	try({ gpar = EMn(data=data, gpar=gpar, label= label, mtol=mtol, mmax= mmax, covtype= covtype, skewness= skewness )$gpar }, TRUE)

	return(gpar)
}
	
trueclassgpar <- function(
	data=NULL, n=10, G=NULL, true.labels=NULL, v=1,  mtol=NULL, mmax=NULL, covtype=NULL, skewness=FALSE, dg=NULL 
) {
  #print("igpar3 1")
  lw = true.class
  z = combinewk(weights=matrix(0,nrow=nrow(data),ncol=G), label=lw)
  gpar = rgpar(data=data, g=G, w=z, d=dg)
  
  for (j in 1:n) try({ gpar = EMgrstep(data=data, gpar=gpar, v=v, label=label, w=z, mtol=mtol, mmax=mmax, covtype=covtype, skewness=skewness )}, TRUE)
  try({ gpar = EMn(data=data, gpar=gpar, label= label, mtol=mtol, mmax= mmax, covtype= covtype, skewness= skewness )$gpar }, TRUE)
  
  return(gpar)
}






	

#logbesselKv <- function(x, y) { log(besselK(x=y, nu=x, expon.scaled=TRUE)) - log(y)}

logbesselKv <- function(x, y) { 
  val = log(besselK(x=y, nu=x, expon.scaled=TRUE)) -y 
  subset = is.nan(val)
  if (any(subset)) val[subset] = besselK.nuAsym(x=y[subset], nu=x[subset], k.max=4, expon.scaled = TRUE, log = TRUE) - y[subset]
  return(subset)
  }

logbesselKv <- function(x, y) { 
 log(besselK(x=y, nu=x, expon.scaled=TRUE)) -y 
}

besselKv    <- function(x, y) { besselK(x=y, nu=x)}

gig <- function(a=NULL,b=NULL,v=NULL) {
	# returns a matrix with dim length(a) x 3
	sab =  sqrt(a*b) 
	kv1 = besselK( sab, nu=v+1, expon.scaled = TRUE)
	kv  =  besselK( sab, nu=v, expon.scaled = TRUE)
	kv12 = kv1/kv

	sb.a = sqrt(b/a)
	w    = kv12*sb.a
	invw = kv12*1/sb.a - 2*v/b
	logw = log(sb.a) + 
	numDeriv::grad( 
		logbesselKv, 
		x=rep(v,length(sab)),
		y=sab, 
		method="Richardson",  
		method.args=list(
			eps=1e-8, d=0.0001, 
			zero.tol=sqrt(.Machine$double.eps/7e-7), 
			r=6, 
			v=2, 
			show.details=FALSE
		)
	)

	val =cbind(w,invw,logw)	

	if (any(is.nan(val))) print( c(sab[170],v[1]))

	return(val)
	}
	
	

gig2 <- function(q1=NULL, q2=NULL, par=NULL) {
  d = length(par$psi)
  p = nrow(par$gam);

  omega1 = par$ol1[1]; lambda1 = par$ol1[2]; 
  omega2 = par$ol2[1]; lambda2 = par$ol2[2]; 
  
  a1 = omega1 + sum( par$beta^2/par$psi )
  b1 = omega1 + q1
  v1 = lambda1- d/2
  
  a2 = omega2 + 0
  b2 = omega2 + q2
  v2 = lambda2-(p-d)/2

  val = list(abc1=gig(b=b1,a=a1, v=v1), abc2=gig(b=b2,a=a2, v=v2))
  return(val)
}



weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )

	



update.maRol <- function(x, par, weights=NULL, alpha.known=NULL, v=1) {
	if (is.null(weights)) weights=rep(1,nrow(x))

	# expectations of w, 1/w, log(w) given x	
	
	### Given Gam we have

	gam.x = x %*% par$gam
	gam.xmu = sweep(gam.x,2,par$mu)
  
    quad1 = apply(sweep( gam.xmu^2, 2, 1/par$psi, "*"),1,sum)
    quad2 = (apply(x^2, 1, sum)  - apply(gam.x^2,1,sum))/par$eta
	
	estep = gig2(quad1, quad2, par=par)
	abc1 = estep$abc1
	abc2 = estep$abc2
	d = length(par$mu); p = nrow(par$gam);

	sumw = sum(weights)
	ABC1 = apply(abc1,2,weighted.sum, wt=weights)/sumw
	ABC2 = apply(abc2,2,weighted.sum, wt=weights)/sumw

	## Now for Gam
    R1 = cov.wt(x, wt=abc1[,2]*weights, center=rep(0,p), method="ML")$cov*ABC1[2]
	R2 = cov.wt(x, wt=abc2[,2]*weights, center=rep(0,p), method="ML")$cov*ABC2[2]/par$eta
	avg.xb = apply(x,2,weighted.sum, wt=abc1[,2]*weights )/sumw
	avg.x  = apply(x,2,weighted.sum, wt=weights )/sumw
	A  = diag(1/par$psi) %*% (outer(par$mu,avg.xb) + outer(par$beta,avg.x) )
	par$gam <- update.gam(par$gam, 1/par$psi, R1, R2, A)

	gam.x = x %*% par$gam
	u1 = (ABC1[2] - abc1[,2])*weights;
	t1 = (ABC1[1]*abc1[,2]-1)*weights; T1 = sum(t1);
	par$mu   = apply(gam.x, 2, weighted.sum, wt=t1)/T1
	par$beta = apply(gam.x, 2, weighted.sum, wt=u1)/T1

	gam.xmu     = sweep(gam.x,2,par$mu)
	a           = apply(gam.xmu^2, 2, weighted.sum, wt=abc1[,2]*weights)/sumw
	avg.gam.xmu = apply(gam.xmu, 2, weighted.sum, wt=weights)/sumw
	par$psi     = a - 2*avg.gam.xmu*par$beta + ABC1[1] * par$beta^2

	quad2   = apply(x^2, 1, sum) - apply(gam.x^2,1,sum) 
	par$eta = sum(quad2*abc2[,2]*weights)/((p-d)*sumw)

	par$ol1 = update.ol(ol= par$ol1, ABC=ABC1, n=1)
	par$ol2 = update.ol(ol= par$ol2, ABC=ABC2, n=1)


	return(par)
	}




update.gam <- function(gam=NULL, ipsi=NULL, R1=NULL, R2=NULL, A=NULL) {
  ## find the minimum for   
  ## tr( ipsi * gam' * R1 * gam )/2 - tr( gam' * R2 * gam )/2  - tr(gam A )
  
  lam1 = eigen(R1, symmetric = TRUE, only.values = TRUE)$values
  G1   = diag(ipsi) %*% t(gam) %*% R1 - max(lam1) * diag(ipsi) %*%  t(gam) 
  
  R2 = -1*R2
  lam2 = eigen(R2, symmetric = TRUE, only.values = TRUE)$values
  G2   = t(gam) %*% R2 - max(lam2)  *  t(gam) 
  
  temp  = -1*(G1 + G2 - A)
  gam = svd(temp)$v %*%  t(svd(temp)$u)
  return(gam)
}





RRlamz <- function(x,lam=NULL,z=1) {
	val =Rlam(x, lam=-lam)+ Rlam(x, lam= lam)
	zval = val - z
	return(zval)
}

Rlam <- function(x, lam=NULL) {
	v1 = besselK(x, nu= lam+1, expon.scaled=TRUE) 
	v0 = besselK(x, nu= lam, expon.scaled=TRUE) 
	val = v1/v0
	return(val)
}


update.ol <- function(ol=NULL, ABC=NULL, n=1) {
	
	for (i in 1:n) {
	if (ABC[3] == 0) {
		ol[2] = 0
	} else {
		bv = numDeriv::grad( logbesselKv, x=ol[2], y=ol[1], method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
		ol[2] = ABC[3]*(ol[2]/bv)
	}

	lam0 = ol[2]
	omg0 = ol[1]
	Rp = Rlam(omg0,lam=+lam0)
	Rn = Rlam(omg0,lam=-lam0)
	f1 = Rp + Rn - (ABC[1]+ABC[2])
	f2 = ( Rp^2 - (2*lam0+1)/omg0 *Rp -1 ) + ( Rn^2 - (2*(-1*lam0)+1)/omg0 *Rn -1 ) 

	if ( ol[1] > f1/f2 ) ol[1] = ol[1] - f1/f2
	}
	return(ol)
}
	
		
	


EMgrstep <- function(data=NULL, gpar=NULL, label=NULL, w=NULL, mtol= NULL, mmax= NULL, covtype=NULL, skewness=TRUE, v=1   ) {
	if (is.null(w)) w = weights(data=data, gpar=gpar,v=v)
	if (!is.null(label)) w = combinewk(weights=w, label= label)

	G= length(gpar$pi);
	d= length(gpar[[1]]$mu);
	for (k in 1:G ) {
		if (skewness) gpar[[k]] = update.maRol(x=data, par=gpar[[k]], weights=w[,k], alpha.known=NULL)
		else gpar[[k]] = update.maRol(x=data, par=gpar[[k]], weights=w[,k], alpha.known=rep(0,d))
	}
	gpar$pi = apply(w,2,mean)

	return(gpar)
	}




	
getall <- function(loglik, eps= epsilon) {
	if (length(loglik) <3) stop("must have at least 3 likelihood values")
	n = length(loglik)
	lm1 = loglik[n]
	lm  = loglik[(n-1)]
	lm_1  = loglik[(n-2)]
	am = (lm1 - lm)/(lm - lm_1)
	lm1.Inf = lm + (lm1 - lm)/(1-am)
	val = lm1.Inf - lm	
	
	continue = TRUE
	if (is.nan(val)) val = 0
	if ( val < eps & val >= 0 )	 continue = FALSE
	return( continue )
	}
	
	
	
EM <- function(data=NULL, gpar0=NULL, G=2, dg = NULL, max.iter=100, epsilon=1e-2, print.loglik=FALSE, label=NULL,  nstart=0, mtol=1e-8, mmax=10, covtype=NULL , skewness= TRUE  ) {
	#Set intrinsic dimension
  if(is.null(dg)) dg = rep(d=floor(0.143*ncol(data)), times=G)
  
  #Set initial parameters
  if (is.null(gpar0)) gpar  = igpar(data=data, g=G, nstart=nstart, mtol= mtol, mmax= mmax, covtype= covtype, skewness= skewness, dg=dg  )
	else gpar = gpar0
	
	#Set number of iterations
	loglik = numeric(max.iter)
	
	#What dat mean?
	for (i in 1:3) {
		loglik[i] = llik(data, gpar)
		gpar = EMgrstep(data=data, gpar=gpar, v=1, label= label, mtol= mtol, mmax= mmax, covtype= covtype, skewness= skewness   )	
	}
	
  #Proceed with EM algorithm to maximize likelihood
	while ( getall(loglik = loglik[1:i], eps=epsilon) & (i < max.iter ) )  {
		i = i+1
 		gpar = EMgrstep(data=data, gpar=gpar, label = label , mtol= mtol, mmax= mmax, covtype= covtype, skewness= skewness  )
		loglik[i] = llik(data, gpar)

	}
	val = list(loglik   = loglik[1:i], 
	           gpar     = gpar, 
	           z        = weights(data=data, gpar= gpar), 
	           map      = MAP(data=data, gpar= gpar, label=label), 
	           skewness = skewness)
	return(val)
	}



EMn <- function(data=NULL, gpar0=NULL, G=2, n=10, label =NULL , nstart=0, mtol=1e-8, mmax=10, covtype= NULL , skewness= TRUE    ) {
	if (is.null(gpar0)) gpar = igpar(data=data, g=G, nstart=nstart, mtol= mtol, mmax= mmax, covtype= covtype, skewness= skewness )
	else gpar  = gpar0
	
	loglik = numeric(n)
	for (i in 1:n) {
 		gpar = EMgrstep(data=data, gpar=gpar, label = label , mtol= mtol, mmax= mmax, covtype= covtype, skewness= skewness   )		
		loglik[i] = llik(data, gpar)
	}
	val = list(loglik   = loglik, 
	           gpar     = gpar, 
	           z        = weights(data=data, gpar= gpar), 
	           map      = MAP(data=data, gpar= gpar, label=label), 
	           skewness = skewness )
	return(val)
	}


update.gam <- function(gam=NULL, ipsi=NULL, R1=NULL, R2=NULL, A=NULL) {
  ## find the minimum for   
  ## tr( ipsi * gam' * R1 * gam )/2 - tr( gam' * R2 * gam )/2  - tr(gam A )
  
  lam1 = eigen(R1, symmetric = TRUE, only.values = TRUE)$values
  G1   = diag(ipsi) %*% t(gam) %*% R1 - max(lam1) * diag(ipsi) %*%  t(gam) 
  
  R2 = -1*R2
  lam2 = eigen(R2, symmetric = TRUE, only.values = TRUE)$values
  G2   = t(gam) %*% R2 - max(lam2)  *  t(gam) 
  
  temp  = -1*(G1 + G2 - A)
  gam = svd(temp)$v %*%  t(svd(temp)$u)
  return(gam)
}  
  
  
  





