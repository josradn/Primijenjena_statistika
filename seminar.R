k_alpha=function(n,alpha)
{
  return ((n-1)*5/(n-5)*qf(1-alpha,df1=5,df2=n-5));
}
b_alpha=function(n,alpha)
{
  return ((qt(1-alpha/10,df=n-1))^2);
}

k1=k_alpha(7,0.05)
b1=b_alpha(7,0.05)
k1
b1
k2=k_alpha(10,0.05)
b2=b_alpha(10,0.05)
k2
b2

#unos podataka
tiroksin=matrix(c(59,85,121,156,191,54,71,90,110,138,56,75,108,151,189,59,85,116,148,177,57,72,97,120,144,52,73,97,116,140,52,70,105,138,171),nrow=7,byrow=T)
tiroksin
tiouracil=matrix(c(61,		 86	,	 109,		 120,		 129,
                   59	,	 80		, 101,		 111	,	 122,
                   53	,	 79	,	 100,		 106,		 133,
                   59	,	 88	,	 100	,	 111	,	 122,
                   51	,	 75		, 101	,	 123	,	 140,
                   51	,	 75		, 92	,	 100	,	 119,
                   56		, 78	,	 95	,	 103	,	 108,
                   58,		 69	,	 93	,	 114	,	 138,
                   46,		 61,		 78	,	 90	,	 107,
                   53	,	 72	,	 89	,	 104	,	 122
),nrow=10,byrow=T)
tiouracil
kontrolna=matrix(c(57,86,114,139,172
                   ,60,93,123,146,177
                   ,52,77,111,144,185
                   ,49,67,100,129,164
                   ,56,81,104,121,151
                   ,46,70,102,131,153
                   ,51,71,94,110,141
                   ,63,91,112,130,154
                   ,49,67,90,112,140
                   ,57,82,110,139,169),nrow=10,byrow=T)
kontrolna

#ocekivanje po vremenima
ocek_tiroksin=c(mean(tiroksin[,1]),mean(tiroksin[,2]),mean(tiroksin[,3]),mean(tiroksin[,4]),mean(tiroksin[,5]))
ocek_tiroksin
ocek_tiouracil=c(mean(tiouracil[,1]),mean(tiouracil[,2]),mean(tiouracil[,3]),mean(tiouracil[,4]),mean(tiouracil[,5]))
ocek_tiouracil
ocek_kontrolna=c(mean(kontrolna[,1]),mean(kontrolna[,2]),mean(kontrolna[,3]),mean(kontrolna[,4]),mean(kontrolna[,5]))
ocek_kontrolna

#kovarijacijska matrica
kovarijacijska_tiroksin=cov(tiroksin, y=tiroksin, use="all.obs")
kovarijacijska_tiroksin
kovarijacijska_tiouracil=cov(tiouracil, y=tiouracil, use="all.obs")
kovarijacijska_tiouracil
kovarijacijska_kontrolna=cov(kontrolna, y=kontrolna, use="all.obs")
kovarijacijska_kontrolna

#pouzdani intervali
A=matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1),nrow=5,byrow=T)
A

#alpha=0.05
n_tiroksin=length(tiroksin[,1])
donje_granice_tiroksin5=c(0,0,0,0,0)
gornje_granice_tiroksin5=c(0,0,0,0,0)
alpha=0.05
for(i in 1:5)
{
  donje_granice_tiroksin5[i]=ocek_tiroksin[i]-sqrt(b_alpha(n_tiroksin,alpha)/n_tiroksin*t(A[,i])%*%kovarijacijska_tiroksin%*%A[,i])
  gornje_granice_tiroksin5[i]=ocek_tiroksin[i]+sqrt(b_alpha(n_tiroksin,alpha)/n_tiroksin*t(A[,i])%*%kovarijacijska_tiroksin%*%A[,i])
}
granice_tiroksin5=matrix(c(donje_granice_tiroksin5,gornje_granice_tiroksin5),nrow=2,byrow=T)
granice_tiroksin5

n_tiouracil=length(tiouracil[,1])
donje_granice_tiouracil5=c(0,0,0,0,0)
gornje_granice_tiouracil5=c(0,0,0,0,0)
alpha=0.05
for(i in 1:5)
{
  donje_granice_tiouracil5[i]=ocek_tiouracil[i]-sqrt(b_alpha(n_tiouracil,alpha)/n_tiouracil*t(A[,i])%*%kovarijacijska_tiouracil%*%A[,i])
  gornje_granice_tiouracil5[i]=ocek_tiouracil[i]+sqrt(b_alpha(n_tiouracil,alpha)/n_tiouracil*t(A[,i])%*%kovarijacijska_tiouracil%*%A[,i])
}
granice_tiouracil5=matrix(c(donje_granice_tiouracil5,gornje_granice_tiouracil5),nrow=2,byrow=T)
granice_tiouracil5

n_kontrolna=length(kontrolna[,1])
donje_granice_kontrolna5=c(0,0,0,0,0)
gornje_granice_kontrolna5=c(0,0,0,0,0)
alpha=0.05
for(i in 1:5)
{
  donje_granice_kontrolna5[i]=ocek_kontrolna[i]-sqrt(b_alpha(n_kontrolna,alpha)/n_kontrolna*t(A[,i])%*%kovarijacijska_kontrolna%*%A[,i])
  gornje_granice_kontrolna5[i]=ocek_kontrolna[i]+sqrt(b_alpha(n_kontrolna,alpha)/n_kontrolna*t(A[,i])%*%kovarijacijska_kontrolna%*%A[,i])
}
granice_kontrolna5=matrix(c(donje_granice_kontrolna5,gornje_granice_kontrolna5),nrow=2,byrow=T)
granice_kontrolna5

#alpha=0.1
donje_granice_tiroksin10=c(0,0,0,0,0)
gornje_granice_tiroksin10=c(0,0,0,0,0)
alpha=0.1
for(i in 1:5)
{
  donje_granice_tiroksin10[i]=ocek_tiroksin[i]-sqrt(b_alpha(n_tiroksin,alpha)/n_tiroksin*t(A[,i])%*%kovarijacijska_tiroksin%*%A[,i])
  gornje_granice_tiroksin10[i]=ocek_tiroksin[i]+sqrt(b_alpha(n_tiroksin,alpha)/n_tiroksin*t(A[,i])%*%kovarijacijska_tiroksin%*%A[,i])
}
granice_tiroksin10=matrix(c(donje_granice_tiroksin10,gornje_granice_tiroksin10),nrow=2,byrow=T)
granice_tiroksin10


donje_granice_tiouracil10=c(0,0,0,0,0)
gornje_granice_tiouracil10=c(0,0,0,0,0)
for(i in 1:5)
{
  donje_granice_tiouracil10[i]=ocek_tiouracil[i]-sqrt(b_alpha(n_tiouracil,alpha)/n_tiouracil*t(A[,i])%*%kovarijacijska_tiouracil%*%A[,i])
  gornje_granice_tiouracil10[i]=ocek_tiouracil[i]+sqrt(b_alpha(n_tiouracil,alpha)/n_tiouracil*t(A[,i])%*%kovarijacijska_tiouracil%*%A[,i])
}
granice_tiouracil10=matrix(c(donje_granice_tiouracil10,gornje_granice_tiouracil10),nrow=2,byrow=T)
granice_tiouracil10


donje_granice_kontrolna10=c(0,0,0,0,0)
gornje_granice_kontrolna10=c(0,0,0,0,0)
for(i in 1:5)
{
  donje_granice_kontrolna10[i]=ocek_kontrolna[i]-sqrt(b_alpha(n_kontrolna,alpha)/n_kontrolna*t(A[,i])%*%kovarijacijska_kontrolna%*%A[,i])
  gornje_granice_kontrolna10[i]=ocek_kontrolna[i]+sqrt(b_alpha(n_kontrolna,alpha)/n_kontrolna*t(A[,i])%*%kovarijacijska_kontrolna%*%A[,i])
}
granice_kontrolna10=matrix(c(donje_granice_kontrolna10,gornje_granice_kontrolna10),nrow=2,byrow=T)
granice_kontrolna10

#alpha=0.01
donje_granice_tiroksin1=c(0,0,0,0,0)
gornje_granice_tiroksin1=c(0,0,0,0,0)
alpha=0.01
for(i in 1:5)
{
  donje_granice_tiroksin1[i]=ocek_tiroksin[i]-sqrt(b_alpha(n_tiroksin,alpha)/n_tiroksin*t(A[,i])%*%kovarijacijska_tiroksin%*%A[,i])
  gornje_granice_tiroksin1[i]=ocek_tiroksin[i]+sqrt(b_alpha(n_tiroksin,alpha)/n_tiroksin*t(A[,i])%*%kovarijacijska_tiroksin%*%A[,i])
}
granice_tiroksin1=matrix(c(donje_granice_tiroksin1,gornje_granice_tiroksin1),nrow=2,byrow=T)
granice_tiroksin1


donje_granice_tiouracil1=c(0,0,0,0,0)
gornje_granice_tiouracil1=c(0,0,0,0,0)
for(i in 1:5)
{
  donje_granice_tiouracil1[i]=ocek_tiouracil[i]-sqrt(b_alpha(n_tiouracil,alpha)/n_tiouracil*t(A[,i])%*%kovarijacijska_tiouracil%*%A[,i])
  gornje_granice_tiouracil1[i]=ocek_tiouracil[i]+sqrt(b_alpha(n_tiouracil,alpha)/n_tiouracil*t(A[,i])%*%kovarijacijska_tiouracil%*%A[,i])
}
granice_tiouracil1=matrix(c(donje_granice_tiouracil1,gornje_granice_tiouracil1),nrow=2,byrow=T)
granice_tiouracil1


donje_granice_kontrolna1=c(0,0,0,0,0)
gornje_granice_kontrolna1=c(0,0,0,0,0)
for(i in 1:5)
{
  donje_granice_kontrolna1[i]=ocek_kontrolna[i]-sqrt(b_alpha(n_kontrolna,alpha)/n_kontrolna*t(A[,i])%*%kovarijacijska_kontrolna%*%A[,i])
  gornje_granice_kontrolna1[i]=ocek_kontrolna[i]+sqrt(b_alpha(n_kontrolna,alpha)/n_kontrolna*t(A[,i])%*%kovarijacijska_kontrolna%*%A[,i])
}
granice_kontrolna1=matrix(c(donje_granice_kontrolna1,gornje_granice_kontrolna1),nrow=2,byrow=T)
granice_kontrolna1


#graf
par(mfrow=c(1,3))
time=0:4

#alpha=0.05
plot(time,ocek_tiroksin, ylim=c(50,200),xlab="vrijeme",ylab="masa",main="Tiroksin (0.05)")
for(i in 1:5)
{
  segments(i-1,granice_tiroksin5[1,i],x1=i-1,y1=granice_tiroksin5[2,i])
  
}
plot(time, ocek_tiouracil,ylim=c(40,140),xlab="vrijeme",ylab="masa",main="Tiouracil (0.05)")
for(i in 1:5)
{
  segments(i-1,granice_tiouracil5[1,i],x1=i-1,y1=granice_tiouracil5[2,i])
  
}
plot(time, ocek_kontrolna,ylim=c(40,180),xlab="vrijeme",ylab="masa",main="Kontrolna (0.05)")
for(i in 1:5)
{
  segments(i-1,granice_kontrolna5[1,i],x1=i-1,y1=granice_kontrolna5[2,i])
  
}

#alpha=0.1
plot(time,ocek_tiroksin, ylim=c(50,200),xlab="vrijeme",ylab="masa",main="Tiroksin (0.1)")
for(i in 1:5)
{
  segments(i-1,granice_tiroksin10[1,i],x1=i-1,y1=granice_tiroksin10[2,i])
  
}
plot(time, ocek_tiouracil,ylim=c(40,140),xlab="vrijeme",ylab="masa",main="Tiouracil (0.1)")
for(i in 1:5)
{
  segments(i-1,granice_tiouracil10[1,i],x1=i-1,y1=granice_tiouracil10[2,i])
  
}
plot(time, ocek_kontrolna,ylim=c(40,180),xlab="vrijeme",ylab="masa",main="Kontrolna (0.1)")
for(i in 1:5)
{
  segments(i-1,granice_kontrolna10[1,i],x1=i-1,y1=granice_kontrolna10[2,i])
  
}

#alpha=0.01
plot(time,ocek_tiroksin, ylim=c(40,210),xlab="vrijeme",ylab="masa",main="Tiroksin (0.01)")
for(i in 1:5)
{
  segments(i-1,granice_tiroksin1[1,i],x1=i-1,y1=granice_tiroksin1[2,i])
  
}
plot(time, ocek_tiouracil,ylim=c(40,140),xlab="vrijeme",ylab="masa",main="Tiouracil (0.01)")
for(i in 1:5)
{
  segments(i-1,granice_tiouracil1[1,i],x1=i-1,y1=granice_tiouracil1[2,i])
  
}
plot(time, ocek_kontrolna,ylim=c(40,200),xlab="vrijeme",ylab="masa",main="Kontrolna (0.01)")
for(i in 1:5)
{
  segments(i-1,granice_kontrolna1[1,i],x1=i-1,y1=granice_kontrolna1[2,i])
  
}


#normalnost podataka
par(mfrow=c(1,5))
for(i in 1:5)
{
  plot(qnorm(((1:n_tiroksin)-3/8)/(n_tiroksin+1/4)),sort(tiroksin[,i]),xlab="qnorm",ylab="tiroksin")
  abline(mean(tiroksin[,i]),sd(tiroksin[,i]))
}

for(i in 1:5)
{
  plot(qnorm(((1:n_tiouracil)-3/8)/(n_tiouracil+1/4)),sort(tiouracil[,i]),xlab="qnorm",ylab="tiouracil")
  abline(mean(tiouracil[,i]),sd(tiouracil[,i]))
}

for(i in 1:5)
{
  plot(qnorm(((1:n_kontrolna)-3/8)/(n_kontrolna+1/4)),sort(kontrolna[,i]),xlab="qnorm",ylab="kontrolna")
  abline(mean(kontrolna[,i]),sd(kontrolna[,i]))
}

library("nortest")
lillie.test(tiroksin[,1])
lillie.test(tiroksin[,2])
lillie.test(tiroksin[,3])
lillie.test(tiroksin[,4])
lillie.test(tiroksin[,5])
lillie.test(tiouracil[,1])
lillie.test(tiouracil[,2])
lillie.test(tiouracil[,3])
lillie.test(tiouracil[,4])
lillie.test(tiouracil[,5])
lillie.test(kontrolna[,1])
lillie.test(kontrolna[,2])
lillie.test(kontrolna[,3])
lillie.test(kontrolna[,4])
lillie.test(kontrolna[,5])

#hotellingov test 
n=n_kontrolna
alpha=0.05
c_alpha=k_alpha(n_kontrolna,alpha)
c_alpha
mu0=c(60,80,100,120,140)
ocek_kontrolna
p=5
n
p
T_2=n*t(ocek_kontrolna-mu0)%*%solve(kovarijacijska_kontrolna)%*%(ocek_kontrolna-mu0)
T_2
pvrijednost=1-pf(T_2,df1=p,df2=n-p)
pvrijednost
