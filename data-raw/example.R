n<-1000 # number of observations
m<-2000 # number of genetic markers
set.seed(1)

# covariates
sex<-factor(rbinom(n,1,.5)+1)
activity<-sapply(floor(rnorm(n,10,5)),function(x)max(x,0))
agecategory<-factor(sample(1:5,n,replace=TRUE))

# effects of covariates
b_sex<-c(0,.1)
b_act<- -.2
b_age<-c(0,.15,.15,.3,.5)

# markers - artificial correlated data
pcond<-matrix(c(0.7405,0.0685,0.001,0.069,0.111,0,0.0005,0.0005,0.009),nrow=3,byrow=TRUE)/c(.81,.18,.01)
  # joint genotype distribution for neighbours (giving correlation 0.1) divided by marginal distribution

xg<-matrix(0,nrow=n,ncol=m)
xg[,1]<-c(rep(0,.81*n),rep(1,.18*n),rep(2,.01*n))
for(j in 2:m) xg[,j]<-sapply(1:n,function(i)sample(0:2,1,prob=pcond[xg[i,j-1]+1,1:3]))

# responses
y<-sapply(1:n,function(i)rnorm(1,50+b_sex[sex[i]]+b_act*activity[i]+b_age[agecategory[i]]))
# effect at marker 1:
eff<-c(0,.43,.83)
# y<-y+sapply(1:n,function(i)if(xg[i,1]==2) .1 else 0) # effect for genotype 2 at marker 1
y<-y+eff[xg[,1]+1] # effect for genotype 2 at marker 1

summary(lm(y~sex+activity+agecategory))

# compute test statistics and first and second order correlations
# library(fwerapprox)

result<-scorestatcorr(y~sex+activity+agecategory,xg,2)
result$statistic[1:10]
result$corrs[[1]][1:10]
pvals<-2*pnorm(-abs(result$statistic))
pvals[1:10]

# find alpha_loc controlling FWER at 0.05 level given by order 2 approximation
al<-uniroot(function(a)gamma2(a,result$corrs[[1]])-.95,c(1e-5,5e-4),tol=1e-14)
al2<-al$root
which(pvals<al2)

# al3<-uniroot(function(a)gamma_k(3,a,result$corrs)-.95,c(1e-5,5e-4),tol=1e-14)
# al3$root

# .05/2000 # Bonferroni
# 1-.95^(1/2000) # Sidak

# logistic model
effl<-c(0,.57,2)
yl<-rbinom(n,1,1/(1+exp(-effl[xg[,1]+1])))

resultl<-scorestatcorr(yl~1,xg,2,family=binomial)
resultl$statistic[1:10]
resultl$corrs[[1]][1:10]
pvalsl<-2*pnorm(-abs(resultl$statistic))
pvalsl[1:10]

# find alpha_loc controlling FWER at 0.05 level given by order 2 approximation
all<-uniroot(function(a)gamma2(a,result$corrs[[1]])-.95,c(1e-5,5e-4),tol=1e-14)
al2l<-all$root
al2l
which(pvalsl<al2l)
