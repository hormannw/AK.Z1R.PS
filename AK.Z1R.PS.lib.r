#####################################################################
# code implementing the Algorithms mentioned in the Dingeç & Hörmann Paper 
# the Algos estimate the right tail probabilities of the sum of exchangeable log-normals 
#
# Authors and Copyright owner:  Wolfgang Hörmann and Kemal Dinçer Dingeç

####################################################
qAsKr <- function(zm,s=2,rho=0.5,gam=1.e8,logyn=T){
# q() simulation function of the Asmussen&Kroese Algorithm for the sum of
# exchangable lognormal vectors using the transform with the D-matrix
# zm ... n x (d-1) input-matrix of iid standard normal variates
# s ... standard deviaton of the 
# rho ... correlation
# gam ... threshold value of P(Sumlognormal < gam)
# logyn ... when TRUE the log value of the result is returned
 d<- dim(zm)[2]+1
 n <- dim(zm)[1]
 alfa <- ifelse(rho==0, 1.5-d/2, 
                 (1-sqrt(1-rho*(rho*(d-2)-d+3)) )/rho )
 Dmat.12 <- 1/sqrt(alfa^2+d-2); Dmat.11 <- alfa*Dmat.12

 zm <- Dmat.11*zm + Dmat.12*(matrix(rowSums(zm),nrow=n,ncol=d-1) - zm)# = D(Z)
 z1 <- log(pmax(0,gam-rowSums(exp(s*zm))))/s 
 #else z1 <- (log1p(pmax(-1,-rowSums(exp(s*zm-log(gam)))))+log(gam))/s # code trying to reduce rounding errors for gam > 1.e20
 z2 <-rowMax(zm)
 icv <-inverseCormat(d=d-1,rho=rho) 
 s.cond <- sqrt(1-rho^2*(d-1)*(icv[1]+(d-2)*icv[2])) 
 mu.cond <- rho*(icv[1]+(d-2)*icv[2])*rowSums(zm)
 if(!logyn) return(d*pnorm(-(pmax(z1,z2)-mu.cond)/s.cond))
 return(log(d) + pnorm(-(pmax(z1,z2)-mu.cond)/s.cond,log.p=T))
}# end qAsKr()

########################################################################################################
mult.ortmatMainDirection <- function(zm){
# calculates the product t(ortmat() %*% t(zm))
# that transform the zm matrix with an orthogonal matrix 
# such that the first variable is transformed into direction (1,1,1,...,1)
# fast implementation with O(n.d) operations
#
 uv <- zm[,-1]
 n <- length(uv[,1])
 d <- length(uv[1,])+1
 hv1 <- c( sqrt(1/((1:(d-1))*(2:d))) , 0 )
 hv2 <- c( 0, sqrt((1:(d-1))/(2:d)) )
 civm <- matrix(0,nrow=n,ncol=d);
 civm[,d] <- -hv2[d]*uv[,d-1]
 cusu <- numeric(n)
 for(i in (d-1):2 ){
   cusu <- cusu+hv1[i]*uv[,i]
   civm[,i] <- cusu -hv2[i]*uv[,i-1] 
 }
 civm[,1] <- cusu + hv1[1]*uv[,1]
 civm+zm[,1]/sqrt(d)
}#end mult.ortmat() 

#####################

inverseCormat <- function(d,rho,noMatrix=F){
# returns the inverse of the Correlation matrix with constant correlation rho
# d, rho ... dimension and correlation of the correlation matrix
# noMatrix... TRUE returns only a vector holding the diagonal value and the off-diagonal value
  y <- (1-1/(1-rho))/(1+(d-1)*rho) # off diagonal element
  x <- y +1/(1-rho) # diagonal element
  if(noMatrix) return(c(diag=x,offDiag=y))
  resm <- matrix(y,ncol=d,nrow=d)
  diag(resm) <- x
  resm
}

#####################

rowMax<- function(mat){
# fast implementation of the maxima of all rows of a matrix.
res<- mat[,1]
for(i in 2:dim(mat)[2]){ res <- pmax(res,mat[,i]) }
return(res)
}

#####################

#######################################################################################################
########################################################################################################
########################################################################################################

AK <- function(n=10000,d=4,s=1,rho=0,gam=51){
# called AK_rho in the paper, Asmussen-Kroese estimate using D-Transform
# returns vector with estimate, SE and relErr
# n ... sample size for simulation
# d, s, rho ... parameters of lognormal exchangeable sum
# gam ... threshold value

  zm <- matrix(rnorm((d-1)*n),nrow=n)
  l.qAsKr <- qAsKr(zm=zm,s=s,gam=gam,rho=rho,logyn=T)
  y <- exp(l.qAsKr+log(1.e100)) # +log(1.e100) to avoid underflow
  c(est=est<-mean(y)*1.e-100,SE=SE<-sqrt(var(y)/n)*1.e-100,relErr=SE/est)
}# end AK()
# system.time(print(AK(n=4.e5,d=30,s=1,rho=0.9,gam=1900)))

#####################################################################################################
#####################################################################################################

AK.IS <- function(n=10000,d=4,s=1,rho=0,gam=51){
# uses the estimate AK_rho with standard IS iid Normal with mean shift zstar of Theorem 1 in main direction
# returns vector with estimate, SE and relErr
# n ... sample size for simulation
# d, s, rho ... parameters of lognormal exchangeable sum
# gam ... threshold value

zstar<-function(){
# simple approximately oprimal mean shift used in Theorem 1 of the right tail paper
rho*log(gam)/(s*sqrt(1+rho*(d-2)))
}
  muIS <- zstar()*sqrt(d-1)
  zm <- matrix(rnorm((d-1)*n),nrow=n)
  zm[,1] <- zm[,1]+muIS
  logw <- muIS*(muIS-2*zm[,1])*0.5
  zm <- mult.ortmatMainDirection(zm)
  l.qAsKr <- qAsKr(zm=zm,s=s,gam=gam,rho=rho,logyn=T)
  y <- exp(l.qAsKr+logw+log(1.e100)) # +log(1.e100) to avoid underflow
  c(est=est<-mean(y)*1.e-100,SE=SE<-sqrt(var(y)/n)*1.e-100,relErr=SE/est)
}# end AK.IS()
#system.time(print(AK.IS(n=4.e5,d=30,s=1,rho=0.9,gam=1900)))
# system.time(print(AK(n=4.e5,d=30,s=1,rho=0.9,gam=1900)))

#####################################################################################################
#####################################################################################################
#####################################################################################################
# functions required for AK.Z1R and AK.Z1R.PS


rchiTDR3<- function(uv,df,constrFac=1,chiTDRl=NULL,onlySetupyn=F){
# constructs the table mountain hat (ie, TDR-3 hat) of the chi distribution with df degrees of freedom
# and generates a sample y of that TDR3-hat using the uniform random variates given as vector uv.
# As that sample is used as IS density the function returns a "length(uv) x 3"-matrix holding as 
# first column the sample y, as second column logw (the log of the IS-weights)
#   and as third column the log-values of the table-mountain IS-density 
# uv ... U(0,1) random variates
# df ... degrees of freedom of the chi-distribution
# constrFac = 1 ... construction point is inflection point
# constrFac < 1 ... more weigth in the tails, less in the center
# chiTDRl ... holding the setup information calcuated in this function
#             if it is NULL the setup is executed
#             if it is not empty no set-up is calculate (ie. the value of df is ignored
# onlySetupyn ... FALSE - setup and generation are executed,   T - only set-up is made and the list TDRl is returned 
dchilog <- function(x,df){# log-density of chi distribution
 (df-1)*log(x)-x^2/2 -(df/2-1)*log(2)-lgamma(df/2)
}
ddchilog <- function(x,df){# derivative of log-density of chi distribution
 (df-1)/x - x
}
Wdchi<-function(df){ # inflection points of pdf of chi distribution used as points of contact for TDR3
  sqrt(df+c(-1,1)*sqrt(8*df-7)/2-0.5)
}

if(length(chiTDRl)==0){ # SETUP to generate from the TDR3-density
 m<-sqrt(df-1)
 lm <- dchilog(m,df)
 xlr <- Wdchi(df)
 xl<- m-(m-xlr[1])*constrFac
 xr<- m+ (xlr[2]-m)*constrFac
 lxl <- dchilog(xl,df)
 lxr <- dchilog(xr,df)
 lsl <- ddchilog(xl,df)
 lsr <- ddchilog(xr,df)
 bl <- xl+(lm-lxl)/lsl
 br <- xr+(lm-lxr)/lsr
 Ac <- (br-bl)*exp(lm)
 Ar <- -exp(lxr+lsr*(br-xr))/lsr
 Al <- (exp(lxl+lsl*(bl-xl))-exp(lxl+lsl*(0-xl)))/lsl
 Alcr<- c(Al,Ac,Ar)/(Al+Ac+Ar)
 Atot<- Al+Ac+Ar
 chiTDRl <-list(Alcr=Alcr,Atot=Atot,lm=lm,xl=xl,lxl=lxl,lsl=lsl,bl=bl,xr=xr,lxr=lxr,lsr=lsr,br=br,df=df)
}
if(onlySetupyn) return(chiTDRl)
l<- chiTDRl

y <- uv# generated RV 
# generate left part
Fxl <- exp(-l$lsl*l$bl)
iv <- uv < l$Alcr[1]
y[iv] <- l$bl+log(Fxl+(uv[iv]/l$Alcr[1])*(1-Fxl))/l$lsl
# generate center part
iv <- uv >= l$Alcr[1] & uv < l$Alcr[1]+l$Alcr[2]
y[iv] <- l$bl+(l$br-l$bl)*(uv[iv]-l$Alcr[1])/l$Alcr[2]
# generate right part
iv <- uv >= l$Alcr[1]+l$Alcr[2]
y[iv] <- l$br+log(1-(uv[iv]-l$Alcr[1]-l$Alcr[2])/l$Alcr[3])/l$lsr

logw <- dchilog(y,l$df)-pmin(l$lxl+(y-l$xl)*l$lsl, l$lm, l$lxr+(y-l$xr)*l$lsr )+log(l$Atot)
return(cbind(y=y,logw=logw,g=pmin(l$lxl+(y-l$xl)*l$lsl, l$lm, l$lxr+(y-l$xr)*l$lsr )-log(l$Atot)))
}#rchiTDR3()
#res<- rchiTDR3(uv=(0:99+runif(100))/100,df=9)


#####################################################################################################

simul.condEz1 <- function(nz1=11,z1min=0,z1max=10,n=1.e3,d=10,s=1,rho=0,gam=1.e3,
                             plotyn=F,printyn=F,logyn=T){
# is required for the setup of AK.Z1R.PS
# simulates f_{Z1}(z1) = E(qAsKr()|z1) for z1 = seq(z1min,z1max,length.out=nz1) 
# for the sake of speed it uses common random numbers, 
# it uses also radius IS necessary to obtain a good performance of AK.Z1R.PS !! 
# returns a nz1 x 2 -matrix holding as first column the z1 values evaluated and as
#         second column the density estimate y or its logarithm
# nz1 ... number of z1 values for which f_Z1 is evaluated 
# z1min, z1max ... smallest and largest value of z1
# d, s, rho ... parameters of the lognormal sum
# gam ... threshold value
# plotyn, printyn to turn on plotting and or printing
# logyn ... to select if the log of y or y should be returned


 logqAsKr.zs <- function(zm,zs=0,s=1,gam=10,rho=0){
 # q() simulation function of the Asmussen&Kroese Algorithm for the sum of
 # exchangable lognormal vectors using the transform with the D-matrix
 # this is a specially fast implementation using the same randomness for all z1 (ie z1) values 
 # returns a  "n x length(zs)" matrix of log(q()) values with n the number of rows of zm
 # zm ... n x (d-1) input-matrix of iid standard normal variates
 # zs ... vector of different z1 values 
 # s, rho ... parameters of lognormalSum distribution 
 # gam ... threshold value of P(Sumlognormal < gam)
 
 d<- dim(zm)[2]+1
 n <- dim(zm)[1]
 alfa <- ifelse(rho==0, 1.5-d/2, 
                 (1-sqrt(1-rho*(rho*(d-2)-d+3)) )/rho )
 Dmat.12 <- 1/sqrt(alfa^2+d-2); Dmat.11 <- alfa*Dmat.12
 zm <- Dmat.11*zm + Dmat.12*(matrix(rowSums(zm),nrow=n,ncol=d-1) - zm)# = D(Z)
 zsAdd <- (Dmat.11 + Dmat.12*(d-2))*zs/sqrt(d-1)
 n.zs <- length(zs) 
 rSesz <- rowSums(exp(s*zm))
 z1.zs <- log(pmax(0,gam-(rSesz%o%exp(s*zsAdd))))/s
 z2 <- rowMax(zm)
 z2.zs <- outer(z2, zsAdd,"+")
 icv <-inverseCormat(d=d-1,rho=rho) 
 s.cond <- sqrt(1-rho^2*(d-1)*(icv[1]+(d-2)*icv[2])) 
 rSzm <- rowSums(zm)
 mu.cond.zs <- rho*(icv[1]+(d-2)*icv[2])*outer(rSzm,(d-1)*zsAdd,"+")
 logres.zs <- log(d) + pnorm(-(pmax(z1.zs,z2.zs)-mu.cond.zs)/s.cond,log.p=T) 
 return(logres.zs)
 }# end qAsKr.zs() in simul.condEz1()

z1v<- seq(z1min,z1max,length.out=nz1)
um <- matrix(rnorm((d-2)*n),nrow=n)
um <- um/sqrt(rowSums(um^2))# random directions   
rv.logw <- rchiTDR3(uv=((0:(n-1))+runif(n))/n,df=d-2,constrFac= 0.4)# this IS important for the stable performance of AK.Z1R.PS
zm <- cbind(0,matrix(rv.logw[,1]*um,nrow=n))#zm <- cbind(0,matrix(rv*um,nrow=n)) 
qmlog <- logqAsKr.zs(mult.ortmatMainDirection(zm),zs=z1v,s=s,gam=gam,rho=rho)
qm <- exp(qmlog+rv.logw[,2])    *rep(dnorm(z1v),each=n)#qm <- exp(qmlog)    *rep(dnorm(z1v),each=n)
yv <- colMeans(qm)
y <- rowMeans(qm)
if(plotyn|printyn){
  resv <- c(est=est<-mean(y)*(z1max-z1min), SE = SE <- sqrt(var(y)/n)*(z1max-z1min), relErr=SE/est)
  if(printyn) print(resv)
  if(plotyn){
    dev.new();plot(z1v,log(yv),type="l")
    SEv <- sqrt(apply(qm,2,var)/n)
    lines(z1v,log(yv-2*SEv),lty=3); lines(z1v,log(yv+2*SEv),lty=3)
  }
}
if(logyn)  return(data.frame(cbind(z1v,yv=log(yv))))
return(data.frame(cbind(z1v,yv)))
}# end simul.condEz1()
#res<-simul.condEz1(nz1=51,z1min=0,z1max=12,n=1000,d=30,s=.25,rho=0,gam=55,plotyn=T)


#####################################################################################################
z1TDR3<- function(uv,nsetup=1500,TDRl=NULL,onlySetupyn=F,constrRatio=0.5,d=60,s=0.25,rho=0,gam=97,plotyn=T,
			nz1.0=101,nz1.1=201,plotPaperYN=F,leftPlotPaper=NULL,rightPlotPaper=NULL,mainPlotPaper="",lwdPlotPaper=1){
# if TDRl==NULL simul.condEz1() is used with sample size nsetup to evaluate the optimal IS density
# and then constructs a TDR3(ie. table mountain)-IS density and makes the TDR3-setup for that IS density.
# in any case it generates a sample of that density using inversion with input uniform vector uv. 
# returns list with first element the matrix y.logw.g with the columns "y", "logw" (IS-weights) "logg" (log-density values)
# and as second element the list TDRl of constants calculated in setup necessary for the generation
#
# uv ... vector of U(0,1) sample
# nsetup ... total sample size used for simul.condEz1() in the setup; must be a multiple of 3
# TDRl ... list returned by an earlier call to that function; if it is handed over all other input variables are ignored
#      ...  NULL means that the setup should be made
# onlySetupyn ... F - setup and genration are executed,   T - only set-up is made and a list with just list TDRl is returned       
# constrRatio ... (r_z in the paper) approximate ratio of f(m)/ f(xl) and  f(m)/f(xr) used for selecting the construction point
# d,s,rho ... parameters of the lognormal sum
# gam=97,
# plotyn=T,
# nz1.0, nz1.1 ... number of z1 values for estimating the f_{Z1} density
#               nz1.0 for first rough estimation, nz1.1 for the second estimation used to construct the table mountain density
#plotPaperYN=F,leftPlotPaper=NULL,rightPlotPaper=NULL,mainPlotPaper="",lwdPlotPaper=1 ... can be ignored
if(length(TDRl)==0){ # SETUP to generate from the TDR3-density
    alfa <- ifelse(rho==0, 1-(d-1)/2, (1-sqrt(1-rho*(rho*(d-2)-d+3)) )/rho )
    Dmat.12 <- 1/sqrt(alfa^2+d-2); Dmat.11 <- alfa*Dmat.12
    c1 <- (Dmat.11 + (d-2)*Dmat.12)
    LRB <-  sqrt(d)*log(gam/d)/(c1*s)

# first the pre-screaning to identify the main region of the marginal distribution
lm= -600# only neceesary for the "try()" to know if it was successful
try( { #start try  necessary to avoid NA for extremely small probabilities!!
   res<-simul.condEz1(nz1=nz1.0,z1min=-30,z1max=max(30,LRB+3),n=nsetup/3,d=d,s=s,rho=rho,gam=gam,plotyn=F,logyn=T)
   iv<- which(res$yv-max(res$yv,na.rm=T)>log(1.e-16))
   z1LB<-res$z1v[max(1,min(iv)-2)]
   z1UB<-res$z1v[min(length(res$z1v),max(iv)+2)]
   if(!is.finite(z1LB*z1UB)){return(0)}
 #the main exploration to find mode and left and right construction points
   res2<-simul.condEz1(nz1=nz1.1,z1min=z1LB,z1max=z1UB,n=nsetup*2/3,d=d,s=s,rho=rho,gam=gam,plotyn=F,logyn=T)
   res <- res2[res2$yv >= -700,]
   step<- (max(res$z1v)-min(res$z1v))/(nz1.1-1)
   iv<-  which(res$yv-max(res$yv)>log(constrRatio)) # used below to find xl and xr
   imax<- which.max(res$yv)
   m <- res$z1v[imax]
   lm <- max(res$yv)                 
 }   ,silent=T ) # end try
 if(lm< -500){return(0)}#as exp(-500) = 7.124576e-218 we can return 0 as we do not try to estimate probs < 1e-200
 il <- min(imax-1,min(iv))-1# il ... first index value for the left construction point 
 repeat({# find a left construction point such that tangent line is really above res$yv (logdiff of 0.1 is tolerated)
         # otherwise try the next point to the left as construction point xl
   xl <- res$z1v[il]
   lxl <- res$yv[il]
   lsl <- (res$yv[il+1]-res$yv[il])/step
   if( sum(lxl+lsl*(res$z1v[1:imax]-xl)+ 0.1 < res$yv[1:imax] ) == 0 ) break 
   else il <- il-1  
   if(il==0)    break # 
 })#end repeat
 if(il>0){ bl <- min(m,xl+(lm-lxl)/lsl) 
    # xl is the left construction point; bl the border between left tail and constant center
 }else{# case left tail log-convex!! 
   # an alternative hat is constucted
   # on (leftmost,xl) it uses a secant between xl and the leftmost point that has function value about f(m)*1.e-16
   # center hat is constant (thus in xl the hat is not continuous!!)
   il <-min(imax-2,min(iv))
   bl <- xl <- res$z1v[il]
   lxl <- res$yv[il]
   lsl <- (res$yv[il]-res$yv[1])/(res$z1v[il]-res$z1v[1])
 }
# start of right side: 
 ir <- max(imax+1,max(iv)) # ir ... first index of right construction point
repeat({ # right side: check if the density is nor more than 0.2 above the hat 
         # otherwise try the next point to the right as construction point xr
   xr <- res$z1v[ir] 
   lxr <- res$yv[ir]
   lsr <- (res$yv[ir+1]-res$yv[ir])/step
   if( sum(lxr+lsr*(res$z1v[-(1:(imax-1))]-xr) +0.2<  res$yv[-(1:(imax-1))]) == 0) break 
   else ir <- ir+1  
   if(ir>length(res$yv)){print("Error z1TDR3()! right tail of marginal density is not log-concave!!!");stop()}  
 })#end repeat                 note that we never experienced the case that this error occured
 # setup to generate random variates of the table mountain function using xl, xr selected above
 br <- max(m,xr+(lm-lxr)/lsr)
 Ac <- (br-bl)*exp(lm)
 Ar <- -exp(lxr+lsr*(br-xr))/lsr
 Al <- (exp(lxl+lsl*(bl-xl)))/lsl
 Atot <- (Al+Ac+Ar)
 Alcr<- c(Al,Ac,Ar)/Atot
 if(plotyn&!plotPaperYN){
  print(paste("bl=",bl,"m=",m,"br=",br))
  print(paste("xl=",xl,"bl=",bl,"m=",m,"br=",br,"xr=",xr))
  print(Alcr)
  zv <- seq(min(res$z1v),max(res$z1v),length.out=500)
  dev.new();plot(zv,ifelse(zv<bl,lxl+(zv-xl)*lsl, ifelse(zv<br,lm, lxr+(zv-xr)*lsr )),type="l",col=2)#,ylim=c(lm-16,lm))
  lines(c(bl,bl),c(-1000,1000),col=2)
  lines(c(m,m),c(-1000,1000),col=3)
  lines(c(br,br),c(-1000,1000),col=4)
  lines(rep(LRB,2),c(-1.e15,1.e15))
  lines(res$z1v,res$yv)
#  lines(zv,-0.5/1*((zv-m))^2+lm,col=3)# asymptotic case for rho=0, variance 1
#  lines(zv,-0.5/(1-0.25)*((zv-m))^2+lm,col=4) # asympto case rho=0.25, variance 0.75
#  lines(zv,-0.5/0.5*((zv-m))^2+lm,col=5) # asympto case rho=0.5
#  lines(zv,-0.5/(1-0.8)*((zv-m))^2+lm,col=6) # case rho=0.8, asympt. variance 0.2
  zv1 <- seq(xl-.1,xr+.1,length.out=500)
  dev.new();plot(zv1,ifelse(zv1<bl,lxl+(zv1-xl)*lsl, ifelse(zv1<br,lm, lxr+(zv1-xr)*lsr )),type="l",col=2)
  lines(res$z1v,res$yv)
 }
 if(plotPaperYN){
  iv<- which(res$z1v > xl-leftPlotPaper & res$z1v < xr+rightPlotPaper)
  plot(res$z1v[iv],res$yv[iv],type="l",xlab="z1",ylab="",main=mainPlotPaper,lwd=lwdPlotPaper)
  zv1 <- seq(min(res$z1v[iv]),max(res$z1v[iv]),length.out=2000)
  lines(zv1,ifelse(zv1<bl,lxl+(zv1-xl)*lsl, ifelse(zv1<br,lm, lxr+(zv1-xr)*lsr )),col=2,lwd=lwdPlotPaper)
 }
 TDRl <-list(Alcr=Alcr,Atot=Atot,lm=lm,xl=xl,lxl=lxl,lsl=lsl,bl=bl,xr=xr,lxr=lxr,lsr=lsr,br=br)

} # end SETUP
if(onlySetupyn) return(list(TDRl=TDRl))
l<- TDRl
# here the generation of the random variates from the table-mountain distribution starts:
y<- uv # y will be the generated values
# generate from part left of the uniform center
iv <- uv < l$Alcr[1]
y[iv] <- l$bl+log((uv[iv]/l$Alcr[1]))/l$lsl
# generate from uniform center part
iv <- uv >= l$Alcr[1] & uv < l$Alcr[1]+l$Alcr[2]
y[iv] <- l$bl+(l$br-l$bl)*(uv[iv]-l$Alcr[1])/l$Alcr[2]
# generate from part right of the uniform center
iv <- uv >= l$Alcr[1]+l$Alcr[2]
y[iv] <- l$br+log(1-(uv[iv]-l$Alcr[1]-l$Alcr[2])/l$Alcr[3])/l$lsr
# y are now the variates of the table mountain distribution
# now logw and logg are calculated
g <- ifelse(y < l$bl, l$lxl+(y-l$xl)*l$lsl, ifelse(y<l$br,l$lm, l$lxr+(y-l$xr)*l$lsr ))
logw <- dnorm(y,log=T)-g+log(l$Atot)
return(list(y.logw.g=cbind(y=y,logw=logw,logg=g),TDRl=TDRl))
 ####
}# end z1TDR3
#
#z1TDR3(uv,nsetup=3750,TDRl=NULL,onlySetupyn=T,d=60,s=0.25,rho=0,gam=97,plotyn=T,constrRatio=1.e-3,nz1.0=101,nz1.1=201)
#z1TDR3(uv,nsetup=3750,TDRl=NULL,onlySetupyn=T,d=20,s=0.25,rho=0.9,gam=119,plotyn=T,constrRatio=0.5,nz1.0=101,nz1.1=201)


###############################################################################
###############################################################################
###############################################################################
AK.Z1R <- function(n=10000,d=4,s=1,rho=0,gam=51,nTDRsetup=3000,constrRatio=0.5,R.TDRFac=0.5,plotyn=F){
# AK IS with TDR3 density (ie. table mountain density) for optimal IS density of Z1 and for chi-distribution for Radius IS
# returns vector with estimate, SE and relErr
# n ... sample size for simulation
# d>=4, s, rho ... parameters of lognormal exchangeable sum
# gam ... threshold value
# nTDRsetup ... number of normal vectors used for set-up simulation; if execution time is very important also 750 or 900 ok
#               must be a multiple of 3, 1/3 is used for the presample to define the domain borders, 
#               2/3 to estimate the marginal density f_Z1 in that domain
# for the radius IS with TDR3-like IS density for the chi-distribution with construction points
# in the inflection points of the density is used.
# R.TDRFac ...    (called c_r in the paper) (default 1) the construction points can be selected closer or 
#                  farer away than the inflection points  R.TDRFac ... 0 ... no radius IS is used
#
# for Z1-IS a TDR-3 hat for the with nTDRsetup common random numbers simulated marginal optimal IS pseudo-dnesity of Z1 is used.
# constrRatio ... (r_z in the paper) =f(xrl)/f(m) (xrl left or right point of contact, should be selected eg. 
#                  0.5 is a good value in our experience
#                     values between 0.2 and 0.8 give similar results                 
# plotyn ... TRUE to plot the marginal density f_Z1 and its upper bound constructed in the set-up

print("hallo")

  um <- matrix(rnorm((d-2)*n),nrow=n)
  um <- um/sqrt(rowSums(um^2)) #... n x (d-2) matrix with n random directions as row vectors
  u1 <- runif(n)
  z1v.logw <-z1TDR3(uv=u1,nsetup=nTDRsetup,TDRl=NULL,d=d,s=s,rho=rho,gam=gam,
                    constrRatio=constrRatio,plotyn=plotyn)
  if(length(z1v.logw)==1 ){ return(c(est=0,SE=0,relERR=1)) }#probability underflow in z1TDR3 
  z1v <- z1v.logw[[1]][,1] 
  u2 <- runif(n)
  if(R.TDRFac > 0){ 
    rv.logw <- rchiTDR3(uv=u2,df=d-2,constrFac=R.TDRFac)
  }else{ # no radius IS
    rv.logw<- cbind(sqrt(qchisq(u2,df=d-2)),0)
  }
  rv <- rv.logw[,1] #
  zm <- mult.ortmatMainDirection(cbind(z1v,um*rv))
  l.qAsKr <- qAsKr(zm=zm,s=s,gam=gam,rho=rho,logyn=T)
  y <- exp(l.qAsKr+z1v.logw[[1]][,2]+rv.logw[,2]+log(1.e100)) # +log(1.e100) to avoid underflow
  c(est=est<-mean(y)*1.e-100,SE=SE<-sqrt(var(y)/n)*1.e-100,relErr=SE/est)
}# end AK.Z1R()
#system.time(print(AK.Z1R(n=4.e5,d=30,s=1,rho=0.9,gam=1900,nTDRsetup=750,R.TDRFac=1,constrRatio=0.5)))


###############################################################################
###############################################################################
###############################################################################
AK.Z1R.PS <- function(nz1=20,nr=nz1,nstr=50,d=4,s=1,rho=0,gam=51,nTDRsetup=3000,constrRatio=0.5,R.TDRFac=0.5,plotyn=F){
# AK_rho with IS with TDR3 densiy for optimal IS density of Z1 and chi-distribution for R and proportional stratification
# same as AK.Z1R but Proportional Stratification for Z1 and R is added
# returns vector with estimate, SE and relErr
# nz1 ... number of different intervals for z1 (ie. main direction) stratification 
# nr ... numer of intervals for radius stratification
# nstr ... number of independent (z1,r) points generated in each of the strata
# uses a total of nstr*nz1*nr evaluations of q()
# d>=4, s, rho ... parameters of lognormal exchangeable sum
# gam ... threshold value
# nTDRsetup ... number of normal vectors used for set-up simulation; if execution time is very important also 750 or 900 ok
#               must be a multiple of 3, 1/3 is used for the presample to define the domain borders, 
#               2/3 to estimate the marginal density f_Z1 in that domain
# for the radius IS with TDR3-like IS density for the chi-distribution with construction points
# in the inflection points of the density is used.
# R.TDRFac ...    (called c_r in the paper) (default 1) the construction points can be selected closer or 
#                  farer away than the inflection points  R.TDRFac ... 
#               0 ... no radius IS is used, not recommended unless gam is in the asymptotic regime!!
# for Z1-IS a TDR-3 hat for the with nTDRsetup common random numbers simulated marginal optimal IS pseudo-dnesity of Z1 is used.
# constrRatio ... (r_z in the paper) =f(xrl)/f(m) (xrl left or right point of contact, should be selected eg. 
#                  0.5 is a good value in our experience
#                     values between 0.2 and 0.8 give similar results                 
# plotyn ... TRUE to plot the marginal density f_Z1 and its upper bound constructed in the set-up

  TDRl <-z1TDR3(uv=NULL,nsetup=nTDRsetup,TDRl=NULL,onlySetupyn=T,d=d,s=s,rho=rho,gam=gam,
              constrRatio=constrRatio,plotyn=plotyn)# set up that constructs the table-mountain IS density
  if(!is.list(TDRl)==1){ return(c(est=0,SE=0,relERR=1)) }#probability underflow in z1TDR3 
  resStr100<- numeric(nstr)
   if(R.TDRFac > 0){ 
     chiTDRl <- rchiTDR3(uv=u2,df=d-2,chiTDRl=NULL,onlySetupyn=T,constrFac=R.TDRFac)# setup for radius IS density
   }  

# in the below for loop independent repetitions of proportional allocation stratification run with 1 point in each stratum 
  m <- nr*nz1
  lb1 <- rep(0:(nz1-1),each=nr)# lower bounds to generate the uniform vector used for z1 strata
  lb2 <- rep(0:(nr-1),nz1)# lower bounds to generate the uniform vector used for r strata
  for(i in 1:nstr){
    um <- matrix(rnorm((d-2)*m),nrow=m)
    um <- um/sqrt(rowSums(um^2)) #... m x (d-2) matrix holding in its rows m random directions
    u1 <- (lb1+runif(m))/nz1
    z1v.logw <-z1TDR3(uv=u1,TDRl=TDRl[[1]])
    z1v <- z1v.logw[[1]][,1]  #a+(b-a)*u1
    u2 <- (lb2+runif(m))/nr
    if(R.TDRFac > 0){ 
      rv.logw <- rchiTDR3(uv=u2,chiTDRl=chiTDRl)
    }else{ # no IS, all log-weights therefore 0 
      rv.logw<- cbind(sqrt(qchisq(u2,df=d-2)),0)
    }
    rv <- rv.logw[,1] 
    zm <- mult.ortmatMainDirection(cbind(z1v,um*rv))# transforms zm in main direction
    l.qAsKr <- qAsKr(zm=zm,s=s,gam=gam,rho=rho,logyn=T)
    yp100 <- exp(l.qAsKr+z1v.logw[[1]][,2]+rv.logw[,2]+log(1.e100)) # +log(1.e100) to avoid underflow when calculating the variance
    resStr100[i] <- mean(yp100)
  }
  c(est=est<-mean(resStr100)*1.e-100,SE=SE<-sqrt(var(resStr100)/nstr)*1.e-100,relErr=SE/est)
}# end AK.Z1R.PS()
#system.time(print(AK.Z1R.PS(nz1=36,nr=36,nstr=31,d=30,s=0.25,rho=0.05,gam=95,nTDRsetup=3000,R.TDRFac=0.3,constrRatio=0.5,plotyn=T)))




# very small performance comparison for a "difficult case"
#system.time(print(AK.Z1R.PS(nz1=36,nstr=31,d=30,s=0.25,rho=0.05,gam=95,nTDRsetup=3000,R.TDRFac=0.3,constrRatio=0.5)))
#system.time(print(    AK(n=4.e5,d=30,s=0.25,rho=0.05,gam=95)))
#system.time(print( AK.IS(n=4.e5,d=30,s=0.25,rho=0.05,gam=95)))
#system.time(print(AK.Z1R(n=4.e5,d=30,s=0.25,rho=0.05,gam=95,nTDRsetup=3000,R.TDRFac=0.3,constrRatio=0.5)))
#system.time(print(AK.Z1R.PS(nz1=100,nstr=40,d=30,s=0.25,rho=0.05,gam=95,nTDRsetup=3000,R.TDRFac=0.3,constrRatio=0.5)))


# and here for an asymptotic case
#system.time(print(AK.Z1R.PS(nz1=36,nstr=31,d=30,s=0.25,rho=0.05,gam=295,nTDRsetup=3000,R.TDRFac=0.3,constrRatio=0.5)))
#system.time(print(    AK(n=4.e5,d=30,s=0.25,rho=0.05,gam=295)))
#system.time(print( AK.IS(n=4.e5,d=30,s=0.25,rho=0.05,gam=295)))
#system.time(print(AK.Z1R(n=4.e5,d=30,s=0.25,rho=0.05,gam=295,nTDRsetup=3000,R.TDRFac=0.3,constrRatio=0.5)))
#system.time(print(AK.Z1R.PS(nz1=100,nstr=40,d=30,s=0.25,rho=0.05,gam=295,nTDRsetup=3000,R.TDRFac=0.3,constrRatio=0.5)))

