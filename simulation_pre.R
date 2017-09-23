#compute average loss of CLIME and RCLIME
library(mvtnorm)
library(flare)
library(Matrix)
source("rob.r")
source("datagenerate.r")
#iter is iteration time
rpre=function(iter=100,N=100,nu=3,P=100,Q=5,Rho=0.6,Type="t", Model="band", Mu=seq(0, 0, length.out=P), S=1 ){
  lossrclime2=rep(0, iter)
  lossclime2=rep(0, iter)
  lossrclimef=rep(0, iter)
  lossclimef=rep(0, iter)
  lossrclimei=rep(0, iter)
  lossclimei=rep(0, iter)
  for (k in 1:iter){
    D=data_generator(n=N,d=nu,p=P,q=Q,rho=Rho,type=Type, mu=Mu, model = Model, s=S) #generate a dataset
    covrob=robcov5(D$data)            #robust cov
    #rclime
    r=sugm(covrob$cov, lambda = NULL, nlambda =20, lambda.min.ratio = 0.1, rho = NULL, method = "clime", sym = "or", shrink=NULL, prec = 1e-4, max.ite = 1e4, standardize = FALSE, perturb = TRUE, verbose = TRUE)
    m=sum(r$icov[[1]]*D$sigma)-determinant(r$icov[[1]])$modulus
    Lambda=1
    for (i in 2:20){
      l=sum(r$icov[[i]]*D$sigma)-determinant(r$icov[[i]])$modulus
      #l=sum(r$icov[[i]]*r$)-determinant(r$icov[[i]])$modulus
      if (l<m){
        m=l
        Lambda=i
      }
    }
    prerob=r$icov[[Lambda]]
    #clime
    c=sugm(D$data, lambda = NULL, nlambda =20, lambda.min.ratio = 0.1, rho = NULL, method = "clime", sym = "or", shrink=NULL, prec = 1e-4, max.ite = 1e4, standardize = FALSE, perturb = TRUE, verbose = TRUE)
    m=sum(c$icov[[1]]*D$sigma)-determinant(c$icov[[1]])$modulus
    lambda=1
    for (i in 2:20){
      l=sum(c$icov[[i]]*D$sigma)-determinant(c$icov[[i]])$modulus
      if (l<m){
        m=l
        lambda=i
      }
    }
    preclime=c$icov[[lambda]]
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
           se_i_rclime=sd(lossrclimei), se_i_clime=sd(lossclimei))
  return (val)
}
#mu=0, band
s1=rpre(N=300, P=300, iter = 20)
sl1=rpre(N=300, P=300, Type = "normal",iter = 20)
s2=rpre(N=200, P=200, iter = 20)
sl2=rpre(N=200, P=200, Type = "normal",iter = 20) 
s3=rpre(N=200, P=300, iter = 20)
sl3=rpre(N=200, P=300, Type = "normal",iter = 20) 

#mu=seq(10, 40, length.out=p),band
s12=rpre(P=60,Mu=seq(10, 40, length.out=P))
s22=rpre(P=100,Mu=seq(10, 40, length.out=P))
s32=rpre(P=200,Mu=seq(10, 40, length.out=P))
sl12=rpre(P=60, Type = "normal", Mu=seq(10, 40, length.out=P))
sl22=rpre(P=100, Type = "normal", Mu=seq(10, 40, length.out=P))
sl32=rpre(P=200, Type = "normal", Mu=seq(10, 40, length.out=P))

#mu=0, nonsparse
s13=rpre(P=60, Model = "nonsparse")
s23=rpre(P=100,Model = "nonsparse")
s33=rpre(P=200,Model = "nonsparse")
sl13=rpre(P=60, Type = "normal",Model = "nonsparse")
sl23=rpre(P=100, Type = "normal",Model = "nonsparse")
sl33=rpre(P=200, Type = "normal",Model = "nonsparse")

s40=rpre(P=300)
sl40=rpre(P=300, Type = "normal")
s41=rpre(P=300,Mu=seq(10, 40, length.out=P))
sl42=rpre(P=300, Type = "normal", Mu=seq(10, 40, length.out=P))
s43=rpre(P=300,Model = "nonsparse")
sl43=rpre(P=300, Type = "normal",Model = "nonsparse")



test_data=band_generator(n=100,d=3,p=100,q=100,rho=0.6,type="t", mu=seq(0, 0, length.out=P), model = "band") 
asdf=r$icov[[50]]
covrob=robcov4(test_data$data) 
r=sugm(covrob$cov, lambda = NULL, nlambda =5, lambda.min.ratio = NULL, rho = NULL, method = "clime", sym = "or", shrink=NULL, prec = 1e-4, max.ite = 1e4, standardize = FALSE, perturb = TRUE, verbose = TRUE)

 