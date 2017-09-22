## find solution to sum_{i=1}^{n}{\rho_{k}(x_{i}-\theta)}=0
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

##find robust mean
##input: data is a n*p matrix representing a random sample of size n and dimention d;
##       k is a number representing threshold
##output: mean of data
robmeank<-function(data,k){
  d<-length(data[1,])
  u<-rep(0,d)
  for (i in 1:d) u[i]<-robsol(data[,i],k)
  return (u)
}
 
##method 1, compute robust estimation of covariance matrix
##input: data is a n*p matrix representing a random sample of size n and dimention d;
##       k is a number representing threshold, here we use one k for all data
##output: estimated covariance matrix
robcovk<-function(data,k){
  d<-length(data[1,])
  u<-rep(0,d)
  cov<-matrix(nrow=d,ncol=d)
  for (i in 1:d) u[i]<-robsol(data[,i],k)
  for (i in 1:d){
    for (j in i:d){
      x<-data[,i]*data[,j]
      cov[i,j] = robsol(x,k)-u[i]*u[j]
      cov[j,i] = cov[i,j]
    }
  }
  return (cov)
}

##method 2, compute robust estimation of covariance matrix
##input: data is a n*p matrix representing a random sample of size n and dimention d;
##k is a number representing threshold, here we use k for estimating sample mean, 
##k^2 for estimating second moment
##output: estimated covariance matrix
robcovk2<-function(data,k){
  d<-length(data[1,])
  u<-rep(0,d)
  cov<-matrix(nrow=d,ncol=d)
  for (i in 1:d) u[i]<-robsol(data[,i],k)
  for (i in 1:d){
    for (j in i:d){
      x<-data[,i]*data[,j]
      cov[i,j] = robsol(x,k^2)-u[i]*u[j]
      cov[j,i] = cov[i,j]
    }
  }
  return (cov)
}

##method 3, compute robust estimation of covariance matrix 
##here we use quantile method to obtain two different thresholds when estianting mean and second moment
##input: data is a n*p matrix representing a random sample of size n and dimention d;
##output: cov is estimated covariance matrix
##mu is estimated mean
##k1 and k2 are two thresholds
robcov3<-function(data,q=0.95){
  k1 =(quantile(data, q, names = F)-quantile(data, 1-q, names = F))/2
  k2 =(quantile(data^2, q, names = F)-quantile(data^2, 1-q, names = F))/2
  d<-length(data[1,])
  u<-rep(0,d)
  cov<-matrix(nrow=d,ncol=d)
  for (i in 1:d) u[i]<-robsol(data[,i],k1)
  for (i in 1:d){
    for (j in i:d){
      x=data[,i]*data[,j]
      cov[i,j] = robsol(x,k2)-u[i]*u[j]
      cov[j,i] = cov[i,j]
    }
  }
  val=list(k1=k1, k2=k2, cov=cov, mu=u)
  return (val)
}

##method 4, compute robust estimation of covariance matrix 
##here we use quantile method to obtain p*(p+1) different thresholds when estianting mean and second moment
##input: data is a n*p matrix representing a random sample of size n and dimention d;
##output: cov is estimated covariance matrix
##mu is estimated mean
##k1 and k2 are thresholds
robcov4<-function(data,q=0.95){
  k1 = apply(data, 2, function(x) (quantile(x, q, names = F)-quantile(x, 1-q, names = F))/2)
  d<-length(data[1,])
  u<-rep(0,d)
  cov<-matrix(nrow=d,ncol=d)
  k2<-matrix(nrow=d,ncol=d)
  for (i in 1:d) u[i]<-robsol(data[,i],k1[i])
  for (i in 1:d){
    for (j in i:d){
      x=data[,i]*data[,j]
      k2[i,j]=(quantile(x, q, names = F)-quantile(x, 1-q, names = F))/2
      cov[i,j] = robsol(x,k2[i,j])-u[i]*u[j]
      cov[j,i] = cov[i,j]
    }
  }
  val=list(k1=k1, k2=k2, cov=cov, mu=u)
  return (val)
}

##method 5, compute robust estimation of covariance matrix 
##here we use quantile method to obtain p different thresholds when estianting mean 
#when computing second moment we demean first and then use quantile method to obtain p^2 different thresholds
##input: data is a n*p matrix representing a random sample of size n and dimention d;
##output: cov is estimated covariance matrix
##mu is estimated mean
##k1 and k2 are thresholds
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


##robust covariance matrix estimation using cross validation to determine tuning parameter k
##we tried both sample covrance matrix and robust covariance matrix as target
robcovcross<-function(data){
  n<-length(data[,1])
  n1<-as.integer(n/2)
  k_tuning<-1     #initiate tuning parameter
  risk<-1000      #initiate minirisk
  for (k in 1:40){
    s<-0
    #20 validation data sets
    for (m in 1:20){ 
      datasample<-data[sample(1:n,n),] #resample
      data1<-datasample[1:n1,]         #training data
      data2<-datasample[(n1+1):n,]     #testing data
      cov1<-robcovk2(data1,k)           #robust estimator
      cov2<-robcovk2(data2,k) 
      #cov2<-(t(data2)-apply(data2,2,mean))%*%t(t(data2)-apply(data2,2,mean))/(n-n1) #sample covariance
      s<-s+norm((cov1-cov2),"F")             ##minimize 2-norm
    }
    s<-s/20
    if (s<risk){
      risk<-s
      k_tuning<-k
    }
  }
  value<-list(tuning=k_tuning, cov=robcovk(data, k_tuning))
  return (value)
}



