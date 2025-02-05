library(readxl)
library(PLRModels)

data <- read_excel("Downloads/data_Indonesia.xlsx")
data <- data[order(data$COVal),]
data

x=data$COVal
y=data$AQIVal


gcv1<-function(y,x,m,l)
{
  a<-min(x)+1
  b<-max(x)-1
  k<-seq(a,b,l)
  v<-length(k)
  n<-length(y)
  Gcv<-matrix(nrow=v,ncol=1)
  Mse<-matrix(nrow=v,ncol=1)
  for (j in 1:v)
  {
    w<-matrix(0,ncol=m+1,nrow=n)
    for (i in 1:m)
      w[,i]<-x^(i-1)
    for (i in m+1)
      w[,i]<-trun(x,k[j],m-1)
    wtw<- t(w) %*% w
    z<- MPL(wtw) 
    beta<- z %*% (t(w) %*% y)
    h<- w %*% z %*% t(w)
    mu<-w%*%beta
    MSE<- t(y-mu) %*% (y-mu)/n
    I<-matrix(0,ncol=n,nrow=n)
    for(i in 1: n)
      I[i,i]<-1
    GCV<-(n^2*MSE)/(sum(diag(I-h)))^2
    Gcv[j]<-GCV
    Mse[j]<-MSE
  }
  R<-matrix(c(k,Gcv,Mse),ncol=3)
  sort.R<-R[order(R[,2]),]
  S<-sort.R[1:10,]
  cat("Untuk spline order",m,"dengan 1 titik knot, diperoleh knot optimal=",S[1,1]," dengan GCV minimum=",S[1,2],"dan MSE =",S[1,3])
  cat("\nBerikut 10 nilai GCV terkecil, nilai MSE dan letak titik knotnya:\n")
  cat("====================================\n")
  cat("  No  Ttk knot   GCV     MSE   \n")
  cat("====================================\n")
  S
}

MPL<-function(x,eps=1e-20)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1])%*%t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-
      xsvd$v[,1:length(diago)]%*%diag(1/diago)%*%t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}

model.spline=function(prediktor,respon,m,knots=c(...))
{
  y<-respon
  n<-length(y)
  k<-length(knots)
  w<-matrix(0, ncol=m+k, nrow=n)
  for (i in 1:m)
    w[,i]<-prediktor^(i-1)
  for(i in (m+1):(m+k))
    w[,i]<-trun(prediktor,knots[i-m],m-1)
  wtw<-t(w)%*%w
  Z<-MPL(wtw)
  beta<-Z%*%t(w)%*%y
  yfits<-w%*%beta
  res<-y-yfits
  MSE<-t(y-yfits)%*%(y-yfits)/n
  I<-matrix(0,ncol=n,nrow=n)
  for(i in 1:n)
    I[i,i]<-1
  h<-w%*%MPL(wtw)%*%t(w)
  GCV<-(n^2*MSE)/(sum(diag(I-h)))^2
  q<-seq(min(prediktor),max(prediktor),length=1000)
  u<-matrix(0,ncol=m+k,nrow=1000)
  cat("\n Spline orde",m)
  cat("\n Titik Knots  = c( ",format(knots),")")
  cat("\n Nilai GCV    = ",format(GCV),
      "\n Nilai MSE    = ",format(MSE),"\n")
  cat("\n ******************************************************************")
  cat("\n      Koefisen         Estimasi")
  cat("\n ******************************************************************")
  for(i in 1:(m+k))
    cat("\n     beta[",i-1,"]          ",format(beta[i]))
  cat("\n ******************************************************************")
  par(mfrow=c(1,1))
  z0=cbind(prediktor,respon)
  z1=z0[order(z0[,1]),]
  x1=z1[,1]
  y1=z1[,2]
  w1<-matrix(0, ncol=m+k, nrow=n)
  for (i in 1:m)
    w1[,i]<-x1^(i-1)
  for(i in (m+1):(m+k))
    w1[,i]<-trun(x1,knots[i-m],m-1)
  yfits1<-w1%*%beta
  plot(x1,y1, type="p",xlim=c(min(prediktor),max(prediktor)),ylim=c(0,300),
       xlab="COVal",ylab="AQIVal")
  par(new=T)
  plot(x1,yfits1, type="l",col="red",
       xlim=c(min(prediktor),max(prediktor)),
       ylim=c(0,300),
       xlab="  ",ylab="  ")
}

trun <- function(data,a,power)
{
  data[data<a] <- a
  (data-a)^power
}


gcv1(y,x,2,1)
model.spline(x,y,2,c(4))


gcv1(y,x,3,1)
model.spline(x,y,3,c(7))
