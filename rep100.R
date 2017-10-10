library(flare)
library(mvtnorm)
library(Matrix)
robsol=function(x, k){
  a<-min(x)
  b<-max(x)
  f = function(u){
    s=ifelse(abs(x-u)<k, x-u, sign(x-u)*k)
    s=sum(s)
    return(s)
  }
  res = uniroot(f, c(a,b))
  return (res$root)
}
robcov5<-function(data, q=0.95){
  d = length(data[1,])
  u = rep(0,d)
  k1 = rep(0,d)
  cov = matrix(nrow=d,ncol=d)
  k2 = matrix(nrow=d,ncol=d)
  for (i in 1:d) {
    k1[i] = (quantile(data[,i], q, names = F)-quantile(data[,i], 1-q, names = F))/2
    u[i] = robsol(data[,i],k1[i])
  }
  for (i in 1:d){
    for (j in i:d){
      x = (data[,i]-u[i])*(data[,j]-u[j])
      k2[i,j] = (quantile(x, q, names = F)-quantile(x, 1-q, names = F))/2
      cov[i,j] = robsol(x, k2[i,j])
      cov[j,i] = cov[i,j]
    }
  }
  return (list(cov=cov, mu=u, k1=k1, k2=k2))
}
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
  else if (type=="mixt"){   
    Z =  matrix(rep(rbinom(n, 1, prob),p),n,p)
    data0 = (rmvt(n, sigma=(d-2)/d*sigma0, df=d) + rep(mu, each=n))*Z + (rmvt(n, sigma=(d-2)/d*s*diag(p), df=d) + rep(mu, each=n)) *(1-Z)
    sigma0=prob^2*sigma0+(1-prob)^2*(s*diag(p)) 
    omega0=solve(sigma0)
  }
  else if (type=="t") data0=rmvt(n, sigma = (d-2)/d*sigma0, df=d) + rep(mu, each=n)
  else if (type=="normal") data0=rmvnorm(n, sigma = sigma0) + rep(mu, each=n)
  else{
    cat("type must be one of t, nomal, mixt or mixnormal")
    return(NULL)
  }
  val=list(omega=omega0, sigma=sigma0, data=data0)
  return (val)
}
rpre=function(iter=100,N=100,nu=3,P=200,Q=5,Rho=0.6,Type="t", Model="band", Mu=seq(0, 0, length.out=200), ss=50, Prob=0.9, thre=0.88){
  lossrclime2=rep(0, iter)
  lossclime2=rep(0, iter)
  lossrclimef=rep(0, iter)
  lossclimef=rep(0, iter)
  lossrclimei=rep(0, iter)
  lossclimei=rep(0, iter)
  for (k in 1:iter){
    #generate a dataset
    D=data_generator(n=N,p=P,type=Type, d=nu,s=ss, prob=Prob,model = Model,q=Q,rho=Rho,mu=Mu)
    #robust cov
    covrob=robcov5(D$data, thre)
    lambda=seq(0.02,0.4,length.out = 10)
    #rclime
    r=sugm(covrob$cov, lambda = lambda, rho = NULL, method = "clime", sym = "or", shrink=NULL, prec = 1e-4, max.ite = 1e4, standardize = FALSE, perturb = TRUE, verbose = TRUE)
    m=norm(D$omega-r$icov[[1]], "2")
    Lambda=1
    for (i in 2:10){
      l=norm(D$omega-r$icov[[i]], "2")
      if (l<m){
        m=l
        Lambda=i
      }
    }
    prerob=r$icov[[Lambda]]
    #clime
    S = cov(D$data, use = "pairwise.complete.obs")
    c=sugm(S, lambda = lambda, rho = NULL, method = "clime", sym = "or", shrink=NULL, prec = 1e-4, max.ite = 1e4, standardize = FALSE, perturb = TRUE, verbose = TRUE)
    m=norm(D$omega-c$icov[[1]], "2")
    Lambda=1
    for (i in 2:10){
      l= norm(D$omega-c$icov[[i]], "2")
      if (l<m){
        m=l
        Lambda=i
      }
    }
    preclime=c$icov[[Lambda]]
    #
    lossrclime2[k]<-norm((prerob-D$omega),"2")
    lossclime2[k]<-norm((preclime-D$omega),"2")
    lossrclimef[k]<-norm((prerob-D$omega),"F")
    lossclimef[k]<-norm((preclime-D$omega),"F")
    lossrclimei[k]<-norm((prerob-D$omega),"i")
    lossclimei[k]<-norm((preclime-D$omega),"i")
  }
  #output average norm loss and sdtandar error under different norms
  val=list(aveloss_2_rclime=mean(lossrclime2), aveloss_2_clime=mean(lossclime2),
           aveloss_F_rclime=mean(lossrclimef), aveloss_F_clime=mean(lossclimef),
           aveloss_i_rclime=mean(lossrclimei), aveloss_i_clime=mean(lossclimei),
           se_2_rclime=sd(lossrclime2), se_2_clime=sd(lossclime2),
           se_F_rclime=sd(lossrclimef), se_F_clime=sd(lossclimef),
           se_i_rclime=sd(lossrclimei), se_i_clime=sd(lossclimei),
           loss_2_rclime=lossrclime2, loss_2_clime=lossclime2,
           loss_F_rclime=lossrclimef, loss_F_clime=lossclimef,
           loss_i_rclime=lossrclimei, loss_i_clime=lossclimei)
  return (val)
}

s1=rpre(iter=100,N=100,nu=3,P=300,Q=5,Rho=0.6,Type="t", Model="band", Mu=seq(0, 0, length.out=300), ss=50, Prob=0.9)
cat(capture.output(print(s1), file="s1.txt"))

s2=rpre(iter=100,N=100,nu=3,P=300,Q=5,Rho=0.6,Type="mixt", Model="band", Mu=seq(0, 0, length.out=300), ss=50, Prob=0.9)
cat(capture.output(print(s2), file="s2.txt"))

s3=rpre(iter=100,N=100,nu=3,P=300,Q=5,Rho=0.6,Type="normal", Model="band", Mu=seq(0, 0, length.out=300), ss=50, Prob=0.9)
cat(capture.output(print(s3), file="s3.txt"))

s10=rpre(iter=100,N=100,nu=3,P=300,Q=5,Rho=0.6,Type="mixnormal", Model="band", Mu=seq(0, 0, length.out=300), ss=50, Prob=0.9)
cat(capture.output(print(s10), file="s10.txt"))

s4=rpre(iter=100,N=100,nu=3,P=100,Q=5,Rho=0.6,Type="t", Model="band", Mu=seq(0, 0, length.out=100), ss=50, Prob=0.9)
cat(capture.output(print(s4), file="s4.txt"))

s5=rpre(iter=100,N=100,nu=3,P=100,Q=5,Rho=0.6,Type="mixt", Model="band", Mu=seq(0, 0, length.out=100), ss=50, Prob=0.9)
cat(capture.output(print(s5), file="s5.txt"))

s6=rpre(iter=100,N=100,nu=3,P=100,Q=5,Rho=0.6,Type="normal", Model="band", Mu=seq(0, 0, length.out=100), ss=50, Prob=0.9)
cat(capture.output(print(s6), file="s6.txt"))

s11=rpre(iter=100,N=100,nu=3,P=100,Q=5,Rho=0.6,Type="mixnormal", Model="band", Mu=seq(0, 0, length.out=100), ss=50, Prob=0.9)
cat(capture.output(print(s11), file="s11.txt"))

s7=rpre(iter=100,N=100,nu=3,P=200,Q=5,Rho=0.6,Type="t", Model="band", Mu=seq(0, 0, length.out=200), ss=50, Prob=0.9)
cat(capture.output(print(s7), file="s7.txt"))

s8=rpre(iter=100,N=100,nu=3,P=200,Q=5,Rho=0.6,Type="mixt", Model="band", Mu=seq(0, 0, length.out=200), ss=50, Prob=0.9)
cat(capture.output(print(s8), file="s8.txt"))

s9=rpre(iter=100,N=100,nu=3,P=200,Q=5,Rho=0.6,Type="normal", Model="band", Mu=seq(0, 0, length.out=200), ss=50, Prob=0.9)
cat(capture.output(print(s9), file="s9.txt"))
 
s12=rpre(iter=100,N=100,nu=3,P=200,Q=5,Rho=0.6,Type="mixnormal", Model="band", Mu=seq(0, 0, length.out=200), ss=50, Prob=0.9)
cat(capture.output(print(s12), file="s12.txt"))

 
