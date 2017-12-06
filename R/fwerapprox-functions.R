#' Score test statistics in generalized linear models, with estimated
#' correlations
#'
#' Computes score test statistics for testing whether each of a large number of
#' coefficients (typically corresponding to genetic markers) in a GLM with canonical
#' link is zero
#' in presence of a smaller number of covariates (typically environmental
#' covariates) not included in the test, and provides estimates of correlations
#' between the test statistics.
#' @param formula Model \code{\link[stats]{formula}} for the null model, i.e.
#'   including environmental covariates.
#' @param xg Matrix of genetic markers (one column for each marker).
#' @param maxorder Maximal order of which correlations of score test statistics
#'   are estimated.
#' @param family Family of GLM passed to \code{\link[stats]{glm}}.
#' @param both If \code{TRUE}, correlations are given both as a list and as a
#'   matrix (see 'Value').
#' @return A list containing the following components: \item{statistic}{A vector
#'   of score test statistics} \item{corrs}{A list of vectors of correlations of
#'   score test statistics. The first component is a vector of first-order
#'   correlations, i.e. between neighbouring markers, the second is a vector of
#'   second-order correlations, i.e. between markers of distance 2, etc. The
#'   number of components is \code{maxorder}. The list is generated if
#'   \code{maxorder} is greater than 0 and less than the number of columns of
#'   \code{xg} or if \code{both} is \code{TRUE}.} \item{corrmatrix}{Estimated
#'   correlation matrix of the score test statistics. The matrix is generated if
#'   \code{maxorder} is equal to the number of columns of \code{xg} or if
#'   \code{both} is \code{TRUE}. In the latter case, any correlations between
#'   markers of distance greater than \code{maxorder} will be set to 0.
#'   \code{scorestatcorr} will attempt to load the \code{Matrix} package and
#'   return \code{corrmatrix} as a sparse banded matrix.}.
#' @examples
#' # Normal model with three environmental covariates:
#' result <- scorestatcorr(y_normal ~ sex + activity + agecategory, xg, 2)
#' result$statistic[1:10]
#' result$corrs[[1]][1:10]
#' pvals <- 2*pnorm(-abs(result$statistic)) # p-values for two-sided test
#' pvals[1:10]
#'
#' # Logistic model without environmental covariates:
#' resultl <- scorestatcorr(y_logistic ~ 1, xg, 2, family = binomial)
#' resultl$statistic[1:10]
#' resultl$corrs[[1]][1:10]
#' pvalsl <- 2*pnorm(-abs(resultl$statistic))
#' pvalsl[1:10]
#' @export
scorestatcorr<-function(formula,xg,maxorder,family=gaussian,both=FALSE){
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
				if(requireNamespace("Matrix",quietly=TRUE))
					corrmatrix<-
						Matrix::bandSparse(m,k=0:maxorder,diagonals=c(list(rep(1,m)),corrs),symmetric=TRUE) else{
				#else
					corrmatrix<-diag(m)
					for(i in 1:maxorder) for(j in 1:(m-i))
						corrmatrix[j,j+i]<-corrmatrix[j+i,j]<-corrs[[i]][j]
				}
		}
	}
	return(list("statistic"=t,"corrs"=corrs,"corrmatrix"=corrmatrix))
}

#' Second-order product-type approximation of multivariate normal probablity
#'
#' Computes a second-order Glaz--Johnson approximation to a multivariate
#' standard normal probability with a given correlation matrix. Gives a
#' familywise error rate level bound in multiple testing for a given local
#' (per-hypothesis) significance level.
#'
#' @param alphaloc Local significance level (se above).
#' @param corr A vector of first-order correlations (i.e., between \eqn{T_j} and
#'   \code{T_{j+1}} for \eqn{j = 1}, \ldots, \eqn{m - 1}).
#' @param tol If \eqn{||\rho| - 1| \le} \code{tol} for a correlation
#'   \eqn{\rho}, then the corresponding factor in the approximation is set to
#'   one.
#' @return \eqn{P(|T_1| < c, \ldots, |T_m| < c)} is approximated for \eqn{(T_1,
#'   \ldots, T_m)} multivariate standard normal with first-order correlations
#'   given by \code{corr}, where \eqn{c} is \code{qnorm(1 - alphaloc/2)}. If
#'   \eqn{(T_1,\ldots T_m)} is a test statistic vector for \eqn{m} hypotheses, a
#'   local significance level of \code{alphaloc}, i.e. rejection of a null
#'   hypothesis if the \eqn{p}-value is less than \code{alphaloc}, will control
#'   familywise error rate at the \code{1 - gamma2(alphaloc, corr)} level.
#' @examples
#' # Normal model with three environmental covariates:
#' result <- scorestatcorr(y_normal ~ sex + activity + agecategory, xg, 2)
#' pvals <- 2*pnorm(-abs(result$statistic))
#' # Find alpha_loc controlling FWER at 0.05 level given by order 2 approximation:
#' al <- uniroot(function(a) gamma2(a, result$corrs[[1]]) - .95, c(1e-5, 5e-4), tol = 1e-14)$root
#' which(pvals < al)
#' 0.05/2000 # Bonferroni
#' 1 - 0.95^(1/2000) # Sidak
#' al # order 2 FWER approximation
#'
#' # Logistic model without environmental covariates:
#' result_l <- scorestatcorr(y_logistic ~ 1, xg, 2, family = binomial)
#' al_l <- uniroot(function(a) gamma2(a, result_l$corrs[[1]]) - .95, c(1e-5, 5e-4), tol = 1e-14)$root
#' @export
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

#' Order \eqn{k} product-type approximation of multivariate normal probablity
#'
#' Computes an order \eqn{k} Glaz--Johnson approximation to a multivariate
#' standard normal probability with a given correlation matrix. Gives a familywise error rate level
#' bound in multiple testing for a given local (per-hypothesis) significance level.
#'
#' The function is quite slow. For second-order approximation, it is better to use \code{\link{gamma2}}.
#'
#' @param k Order of the approximation.
#' @param alphaloc Local significance level (se above).
#' @param corr Either a list of correlations
#'   in which the
#'   first component is a vector of first-order correlations, the second is a vector
#'   of second-order correlations,
#'   up to a vector of \code{k} \eqn{- 1}-order correlations, or a correlation matrix.
#'   Typically the \code{corrs} or \code{corrmatrix} component given by \code{\link{scorestatcorr}}.
#' @param miwasteps The \code{steps} parameter used for \code{Miwa} (\code{\link[mvtnorm]{algorithms}})
#'   in \code{\link[mvtnorm]{pmvnorm}}.
#' @param genz If \code{TRUE}, the \code{GentzBretz} algorithm (\code{\link[mvtnorm]{algorithms}}) is used
#'   for \code{\link[mvtnorm]{pmvnorm}}, otherwise the \code{Miwa} algorithm is used.
#' @return \eqn{P(|T_1| < c, \ldots, |T_m| <
#'   c)} is approximated for \eqn{(T_1, \ldots, T_m)} multivariate standard normal
#'   with covariances (at least up to order \code{k-1}) given by \code{corr}, where \eqn{c} is
#'   \code{qnorm(1 - alphaloc/2)}. If \eqn{(T_1,\ldots T_m)} is a test statistic vector
#'   for \eqn{m} hypotheses, a local significance level of \code{alphaloc}, i.e. rejection of
#'   a null hypothesis if the \eqn{p}-value is less than \code{alphaloc}, will control
#'   familywise error rate at the \code{1 - gamma_k(alphaloc, corr)} level.
#' @examples
#' # Normal model with three environmental covariates:
#' result <- scorestatcorr(y_normal ~ sex + activity + agecategory, xg[,1:200], 2)
#' al <- uniroot(function(a) gamma2(a, result$corrs[[1]]) - .95, c(1e-5, 5e-4), tol = 1e-14)$root
#' al3 <- uniroot(function(a) gamma_k(3, a, result$corrs) - .95, c(1e-5, 5e-4), tol = 1e-14)$root
#' 0.05/2000 # Bonferroni
#' 1 - 0.95^(1/2000) # Sidak
#' al # order 2 FWER approximation
#' al3 # order 3 FWER approximation
#' @export
gamma_k<-function(k,alphaloc,corr,miwasteps=4096,genz=FALSE){
  if(!requireNamespace("mvtnorm",quietly=TRUE)) stop("Please install mvtnorm package",call=FALSE)
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
				temp<-mvtnorm::pmvnorm(rep(-lim,k),rep(lim,k),sigma=corrmat)[[1]]
			} else
				temp<-mvtnorm::pmvnorm(rep(-lim,k),rep(lim,k),sigma=corrmat,alg=mvtnorm::Miwa(miwasteps))[[1]]
		}
		if(temp>1) temp<-1
		prod<-prod*temp
		if(prod-1>1e-7)cat("Numerator >= 1",i,temp,prod,"\n")
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
						temp<-mvtnorm::pmvnorm(rep(-lim,k-1),rep(lim,k-1),sigma=corrmat)[[1]]
					} else
						temp<-mvtnorm::pmvnorm(rep(-lim,k-1),rep(lim,k-1),sigma=corrmat,alg=mvtnorm::Miwa(miwasteps))[[1]]
				}
			}
			if(temp>1) temp<-1
			prod<-prod/temp
			if(prod-1>1e-7)cat("Warning: Denominator >= 1",i,temp,prod,"\n")
		}
	}
  prod
}
