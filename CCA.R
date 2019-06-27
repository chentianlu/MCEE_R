install.packages("pscl")
library("pscl")
data<-read.table('C:\\Users\\thinker\\Desktop\\Rpackage\\关节炎气相.csv',sep=',',header=TRUE)
kk<-0.1
rr<-0.3
####adjustment parameter
numrow<-ncol(data)
Yfactor<-as.matrix(data[,1])
Ffactor<-as.matrix(data[,2])
mfactor<-as.matrix(data[,4:numrow])
nc<-ncol(mfactor)
#####Find F-related metabolites
Frelated<-cancor(Ffactor,mfactor)
nrow<-ncol(mfactor)
Frelatedcor<-array(1:nrow,dim=c(nrow,2))
Frelatedcor[,2]<-Frelated$ycenter
Frelatedcor<-Frelatedcor[order(Frelatedcor[,2],decreasing=F),]
Fnumber<-floor(rr%*%nrow)
Fmetabolite<-as.matrix(Frelatedcor[1:Fnumber,1])
####Find Y-related metabolites
Yrelated<-cancor(Yfactor,mfactor)
Yrelatedcor<-array(1:nrow,dim=c(nrow,2))
Yrelatedcor[,2]<-Yrelated$ycenter
Yrelatedcor<-Yrelatedcor[order(Yrelatedcor[,2],decreasing=F),]
Ynumber<-floor(kk%*%nrow)
Ymetabolite<-as.matrix(Yrelatedcor[1:Ynumber,1])
###adjusted-metabolites
admetabolite<-as.matrix(setdiff(Fmetabolite,Ymetabolite))
admetabolite2<-sort(admetabolite)
####remove F-related component
adjustmatrix<-mfactor[,admetabolite2]
unadjustnumber<-setdiff(seq(from=1, to=nc),admetabolite2)
unadjustmatrix<-mfactor[,unadjustnumber]
nn<-nrow(admetabolite)
parameter<-array(0,dim=c(nn,4))
pvalue<-array(0,dim=c(nn,4))
for(i in 1:nn){
  ####Guassian
  m<-try(glm(adjustmatrix[,i]~Ffactor+Yfactor,family=gaussian()),silent = TRUE)
  m0<-try(m[["coefficients"]][["Ffactor"]],silent = TRUE)
  m2<-try(summary(m),silent=TRUE)
  m3<-try(m2[["coefficients"]][2,4],silent = TRUE)
  m4<-try(as.matrix(pR2(m))[5,],silent = TRUE)
  pvalue[i,1]<-try(m4,silent = TRUE)
  parameter[i,1]<-try(m0,silent=TRUE)
  ###inverse.gaussian
  j<-try(glm(adjustmatrix[,i]~Ffactor+Yfactor,family=inverse.gaussian()),silent = TRUE)
  j0<-try(j[["coefficients"]][["Ffactor"]],silent = TRUE)
  j2<-try(summary(j),silent=TRUE)
  j3<-try(j2[["coefficients"]][2,4],silent = TRUE)
  j4<-try(as.matrix(pR2(j))[5,],silent = TRUE)
  pvalue[i,2]<-try(j4,silent = TRUE)
  parameter[i,2]<-try(j0,silent=TRUE)
  ###poosson
  k<-try(glm(adjustmatrix[,i]~Ffactor+Yfactor,family=poisson()),silent = TRUE)
  k0<-try(k[["coefficients"]][["Ffactor"]],silent = TRUE)
  k2<-try(summary(k),silent=TRUE)
  k3<-try(k2[["coefficients"]][2,4],silent = TRUE)
  k4<-try(as.matrix(pR2(k))[5,],silent = TRUE)
  pvalue[i,3]<-try(k4,silent = TRUE)
  parameter[i,3]<-try(k0,silent=TRUE)
  ###gamma
  q<-try(glm(adjustmatrix[,i]~Ffactor+Yfactor,family=Gamma()),silent = TRUE)
  q0<-try(q[["coefficients"]][["Ffactor"]],silent = TRUE)
  q2<-try(summary(q),silent=TRUE)
  q3<-try(q2[["coefficients"]][2,4],silent = TRUE)
  q4<-try(as.matrix(pR2(q))[5,],silent = TRUE)
  pvalue[i,4]<-try(q4,silent = TRUE)
  parameter[i,4]<-try(q0,silent=TRUE)
}
bestmethod<-array(0,dim=c(nn,1))
for(i in 1:nn){
  bestmethod[i,]<-which.max(pvalue[i,])
}
bestpara<-array(0,dim=c(nn,1))
for(i in 1:nn){
  rr<-bestmethod[i,]
  bestpara[i,]<-parameter[i,rr]
}
####Remove the F-related component
bestpara<-apply(bestpara,2,as.numeric)
for(i in 1:nn){
  para<-bestpara[i,]
  adjustmatrix[,i]<- adjustmatrix[,i]-Ffactor*para
}
allmatrix<-cbind(Yfactor,adjustmatrix,unadjustmatrix)
