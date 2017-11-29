#' Score test statistics in generalized linear models, with estimated correlations
scorestatcorr<-function(formula,xg,maxorder,family="gaussian",both=FALSE){
	corrs<-corrmatrix<-NA
	rownames(xg)<-colnames(xg)<-NULL
	m<-dim(xg)[2]
	fit<-glm(formula,family=family)
	phi<-summary(fit)$dispersion

	xe<-model.matrix(fit)
	muhate<-fit$fitted.values

	u<-as.vector(crossprod(xg,fit$y-muhate))
	lambda<-phi*fit$family$variance(muhate)

	lambdaxg<-lambda*xg
	veg<-crossprod(xe,lambdaxg)
	veeinvveg<-solve(crossprod(xe,lambda*xe))%*%veg

	sd<-sqrt(colSums(xg*lambdaxg)-colSums(veg*veeinvveg))
	t<-u/sd

	if(both||maxorder<m-1) corrs<-vector("list",min(maxorder,m-1))

	if(maxorder>0){
		if(maxorder>=m-1){
			corrmatrix<-cov2cor(crossprod(xg,lambdaxg)-crossprod(veg,veeinvveg))
			if(both) for(i in 1:(m-1)) corrs[[i]]<-diag(corrmatrix[1:(m-i),(i+1):m,drop=FALSE])
		} else{
			for(i in 1:maxorder)
				corrs[[i]]<-
					(colSums(xg[,1:(m-i),drop=FALSE]*lambdaxg[,(i+1):m,drop=FALSE])-
							colSums(veg[,1:(m-i),drop=FALSE]*veeinvveg[,(i+1):m,drop=FALSE]))/
						(sd[1:(m-i)]*sd[(i+1):m])
			if(both)
				if(require(Matrix,quietly=TRUE))
					corrmatrix<-
						bandSparse(m,k=0:maxorder,diagonals=c(list(rep(1,m)),corrs),symmetric=TRUE) else{
				#else
					corrmatrix<-diag(m)
					for(i in 1:maxorder) for(j in 1:(m-i))
						corrmatrix[j,j+i]<-corrmatrix[j+i,j]<-corrs[[i]][j]
				}
		}
	}
	return(list("statistic"=t,"corrs"=corrs,"corrmatrix"=corrmatrix))
}

gamma2<-function(alphaloc,corr,tol=1e-7){ # fast
  const<-sqrt(2/pi)
  m<-length(corr) # one less than the number of markers
  b<-1-alphaloc
  lim<-qnorm(1-alphaloc/2)
  prod<-1
  for(i in 1:m){
  	corr2<-corr[i]
  	if(abs(abs(corr2)-1)>tol){
      const2<-1/sqrt(1-corr2^2)
      prod<-prod*(1-const/b*integrate(function(x)exp(-x^2/2)*pnorm((corr[i]*x-lim)*const2),-lim,lim)$value)
    }
  }
  b*prod
}

gamma_k<-function(k,alphaloc,corr,miwasteps=4096,genz=FALSE){
	require(mvtnorm,quietly=TRUE)
	if(class(corr)=="list") m<-length(corr[[1]])+1 else m<-dim(corr)[1]
	b<-1-alphaloc
	const<-sqrt(2/pi)
	zeromat<-matrix(0,nrow=k,ncol=k)
	prod<-1
	lim<-qnorm(1-alphaloc/2)
	for(i in k:m){
		sing<-FALSE
		if(class(corr)=="list"){
			corrmat<-zeromat
			for(j in 2:k){
				index<-matrix(c(1:(k-j+1),j:k),ncol=2)
				corrmat[index]<-corr[[j-1]][(i-k+1):(i-j+1)]
			}
			corrmat<-corrmat+t(corrmat)+diag(k)
		} else
			corrmat<-as.matrix(corr[(i-k+1):i,(i-k+1):i])
		if(min(eigen(corrmat,only.values=TRUE)$values)<1e-10){
			temp<-1
			sing<-TRUE
		} else {
			if(genz){
				temp<-pmvnorm(rep(-lim,k),rep(lim,k),sigma=corrmat)[[1]]
			} else
				temp<-pmvnorm(rep(-lim,k),rep(lim,k),sigma=corrmat,alg=Miwa(miwasteps))[[1]]
		}
		if(temp>1) temp<-1
		prod<-prod*temp
		if(prod-1>1e-7)cat("Numerator ≥ 1",i,temp,prod,"\n")
		if(i>k){
			if(sing){
				temp<-1
			} else {
				if(k==3){
					corr2<-corrmat[1,2]
					const2<-1/sqrt(1-corr2^2)
					temp<-b-const*integrate(function(x)exp(-x^2/2)*pnorm((corr2*x-lim)*const2),-lim,lim)$value
				} else {
					corrmat<-corrmat[1:(k-1),1:(k-1)]
					if(genz){
						temp<-pmvnorm(rep(-lim,k-1),rep(lim,k-1),sigma=corrmat)[[1]]
					} else
						temp<-pmvnorm(rep(-lim,k-1),rep(lim,k-1),sigma=corrmat,alg=Miwa(miwasteps))[[1]]
				}
			}
			if(temp>1) temp<-1
			prod<-prod/temp
			if(prod-1>1e-7)cat("Warning: Denominator ≥ 1",i,temp,prod,"\n")
		}
	}
  prod
}
