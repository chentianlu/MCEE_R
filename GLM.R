install.packages("pscl")
library("pscl")
data<-read.table('C:\\Users\\thinker\\Desktop\\Rpackage\\关节炎气相.csv',sep=',',header=TRUE)
ss<-0.2####adjustment parameter
numrow<-ncol(data)
Yfactor<-as.matrix(data[,1])
Ffactor<-as.matrix(data[,2])
mfactor<-as.matrix(data[,4:numrow])
nc<-ncol(mfactor)
######The mtabolits affected by confounding factors
pmatrix<-array(0,dim=c(nc,4))
for(i in 1:nc){
  ####Guassian
  a<-try(glm(Ffactor~mfactor[,i],family=gaussian()),silent = TRUE)
  a2<-try(summary(a),silent=TRUE)
  a3<-try(a2[["coefficients"]][2,4],silent = TRUE)
  pmatrix[i,1]<-try(a3,silent = TRUE)
  ###inverse.gaussian
  c<-try(glm(Ffactor~mfactor[,i],family=inverse.gaussian()),silent = TRUE)
  c2<-try(summary(c),silent=TRUE)
  c3<-try(c2[["coefficients"]][2,4],silent = TRUE)
  pmatrix[i,2]<-try(c3,silent = TRUE)
  ###poosson
  d<-try(glm(Ffactor~mfactor[,i],family=poisson()),silent = TRUE)
  d2<-try(summary(d),silent=TRUE)
  d3<-try(d2[["coefficients"]][2,4],silent = TRUE)
  pmatrix[i,3]<-try(d3,silent = TRUE)
  ###logistic
  b<-try(glm(Ffactor~mfactor[,i],family=Gamma()),silent = TRUE)
  b2<-try(summary(b),silent=TRUE)
  b3<-try(b2[["coefficients"]][2,4],silent = TRUE)
  pmatrix[i,4]<-try(b3,silent = TRUE)
}
y=apply(pmatrix,2,as.numeric)
y[is.na(y)] <- 1
pmatrix1<-array(1,dim=c(nc,4))
for(aa in 1:nc){
  for(bb in 1:4){
    if (y[aa,bb]<0.01){
      pmatrix1[aa,bb]<-y[aa,bb]
    }
  }
}
min1<-matrix(apply(pmatrix1,1,min))###每行最小值
min2<-array(1:nc,dim=c(nc,1))
ccc<-array(0,dim=c(nc,1))
for(ii in 1:nc){
  ccc[ii,1]<-which.min(pmatrix1[ii,])
}
min<-cbind(min2,min1,ccc)
min<-min[order(min[,2],decreasing=F),]
min<-subset(min,min[,2]<0.01)
######metabolited affected by Y
pmatrix2<-array(0,dim=c(nc,4))
for(i in 1:nc){
  ####Guassian
  e<-try(glm(mfactor[,i]~Yfactor,family=gaussian()),silent = TRUE)
  e2<-try(summary(e),silent=TRUE)
  e3<-try(e2[["coefficients"]][2,4],silent = TRUE)
  pmatrix2[i,1]<-try(e3,silent = TRUE)
  ###inverse.gaussian
  f<-try(glm(mfactor[,i]~Yfactor,family=inverse.gaussian()),silent = TRUE)
  f2<-try(summary(f),silent=TRUE)
  f3<-try(f2[["coefficients"]][2,4],silent = TRUE)
  pmatrix2[i,2]<-try(f3,silent = TRUE)
  ###poosson
  g<-try(glm(mfactor[,i]~Yfactor,family=poisson()),silent = TRUE)
  g2<-try(summary(g),silent=TRUE)
  g3<-try(g2[["coefficients"]][2,4],silent = TRUE) 
  pmatrix2[i,3]<-try(g3,silent = TRUE)
  ###gamma
  h<-try(glm(mfactor[,i]~Yfactor,family=Gamma()),silent = TRUE)
  h2<-try(summary(h),silent=TRUE)
  h3<-try(h2[["coefficients"]][2,4],silent = TRUE)
  pmatrix2[i,4]<-try(h3,silent = TRUE)
}
y2=apply(pmatrix2,2,as.numeric)
y2[is.na(y2)] <- 1
pmatrix3<-array(1,dim=c(nc,4))
for(aa in 1:nc){
  for(bb in 1:4){
    if (y[aa,bb]<0.01){
      pmatrix3[aa,bb]<-y[aa,bb]
    }
  }
}
min3<-matrix(apply(pmatrix3,1,min))###每行最小值
min4<-array(1:nc,dim=c(nc,1))
ccc1<-array(0,dim=c(nc,1))
for(ii in 1:nc){
  ccc1[ii,1]<-which.min(pmatrix3[ii,])
}
minY<-cbind(min4,min3,ccc1)
minY<-minY[order(minY[,2],decreasing=F),]
minY<-subset(minY,minY[,2]<0.01)
###select metabolites
number<-nrow(minY)*0.33
number<-floor(number)
minY<-minY[1:number,]            
admetabolite<-as.matrix(setdiff(min[,1],minY[,1]))
admetabolite2<-sort(admetabolite)
####adjustment metabolites
adjustmatrix<-mfactor[,admetabolite2]
sunshine<-seq(from=1,to=nc)
unadjustnumber<-setdiff(sunshine,admetabolite2)
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
allmatrix<-cbind(Yfactor,adjustmatrix,unadjustmatrix)####Remove F-related matrix
write.table (allmatrix, file ="C:\\Users\\thinker\\Desktop\\Rpackage", sep ="", row.names =TRUE, col.names =F)