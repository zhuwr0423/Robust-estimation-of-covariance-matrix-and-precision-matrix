#generate data
#input:n=sample size, d=df of t dis, p=dimention, q=band size, rho=band parameter, s=scale of precision matrix
##type="t" or "normal", model="band" or "dense", mu=the vector of noncentrality parameters of length n,
#\sigma_{ij}=\rho^{|i-j|}1{|i-j|<q}
#output:sigma, omega, data
library(mvtnorm)
data_generator=function(n=100,d=3,p=200,q=5,rho=0.6,type="t", model="band", mu=seq(0, 0, length.out=p), s=20){
  omega0=matrix(nrow = p, ncol = p)
  if (model=="band"){
    for (i in 1:p){
      for (j in i:p){
        omega0[i,j]=ifelse((j-i)<q, rho^(j-i), 0)
        omega0[j,i]=omega0[i,j]
      }
    }
  }
  else if (model=="dense"){
    omega0=matrix(rep(0.5, p*p),p,p)
    diag(omega0)=rep(1, p)
  }
  else{
    cat("model must be one of band or dense")
    return(NULL)
    }
  #omega0=omega0*s
  #sigma0=solve(omega0)
  sigma0=cov2cor(solve(omega0))
  omega0=solve(sigma0)
  if (type=="t") data0=rmvt(n, sigma = sigma0, df=d) + rep(mu, each=n)
  if (type=="normal") data0=rmvnorm(n, sigma = sigma0) + rep(mu, each=n)
  val=list(omega=omega0, sigma=sigma0, data=data0)
  return (val)
}
 