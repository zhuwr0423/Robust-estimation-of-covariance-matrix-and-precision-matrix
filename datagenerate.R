#genarate dataset 

#Inpute:
#n:sample size; 
#p:dimension of random vector; 
#mu:the vector of noncentrality parameters of length n,
#type:type of distribution, we have "t", "mixt", "normal", "mixnormal"
#d:df of t distribution; 
#s:strength of outliers in mixture type
#prob: probability of nonoutliers in mixture type
#model: model for precison matrix, we have "band"(\sigma_{ij}=\rho^{|i-j|}1{|i-j|<q});
#"band2"(\sigma_{ij}=\rho^{|i-j|}1{|i-j|<q}*1.05^{j})
#"rand"(diagonal elements are 1, off-diagonal elements are generated independently and equals 0.15 with probability 0.08 or 0 with probability 0.92.)
#"dense"(diagonal elements are 1, off-diagonal elements are 0.5)
#q:banding size; 
#rho:banding strength; 

#Output:omega, sigma, data

library(mvtnorm)
data_generator=function(n=100,p=200,type="t",d=3, s=50, prob=0.9, model="band",q=5,rho=0.6, mu=seq(0, 0, length.out=200)){
  omega0=matrix(nrow = p, ncol = p)
  if (model=="band"){
    for (i in 1:p){
      for (j in i:p){
        omega0[i,j]=ifelse((j-i)<q, rho^(j-i), 0)
        omega0[j,i]=omega0[i,j]
      }
    }
  }
  else if (model=="band2"){
    for (i in 1:p){
      for (j in i:p){
        omega0[i,j]=ifelse((j-i)<q, rho^(j-i)*1.05^(j), 0)
        omega0[j,i]=omega0[i,j]
      }
    }
  }
  else if (model=="rand"){
    omega0 = diag(p)
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        omega0[i,j]=rbinom(1,1,0.08)*0.15
        omega0[j,i]=omega0[i,j]
      }
    }
  }
  else if (model=="dense"){
    omega0=matrix(rep(0.5, p*p),p,p)
    diag(omega0)=rep(1, p)
  }
  else{
    cat("model must be one of band, band2, rand or dense")
    return(NULL)
  }
  sigma0=cov2cor(solve(omega0))
  omega0=solve(sigma0)
  if (type=="mixnormal"){   
    Z =  matrix(rep(rbinom(n, 1, prob),p),n,p)
    data0 = (rmvnorm(n, sigma=sigma0) + rep(mu, each=n))*Z+ (rmvnorm(n, sigma=s*diag(p)) + rep(mu, each=n))*(1-Z)
    sigma0=prob^2*sigma0+(1-prob)^2*(s*diag(p)) 
    omega0=solve(sigma0)
  }
  if (type=="mixt"){   
    Z =  matrix(rep(rbinom(n, 1, prob),p),n,p)
    data0 = (rmvt(n, sigma=(d-2)/d*sigma0, df=d) + rep(mu, each=n))*Z + (rmvt(n, sigma=(d-2)/d*s*diag(p), df=d) + rep(mu, each=n)) *(1-Z)
    sigma0=prob^2*sigma0+(1-prob)^2*(s*diag(p)) 
    omega0=solve(sigma0)
}
if (type=="t") data0=rmvt(n, sigma = (d-2)/d*sigma0, df=d) + rep(mu, each=n)
if (type=="normal") data0=rmvnorm(n, sigma = sigma0) + rep(mu, each=n)
val=list(omega=omega0, sigma=sigma0, data=data0)
return (val)
}
