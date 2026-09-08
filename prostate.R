# TO INSTALL KEGGGRAPH USE THE CODE BELOW
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGGgraph")

##############################################
##############################################
##############################################

library(KEGGgraph)
library(sna)

source("BVSG.r")

pathway=read.csv(file = "pathway.csv",head=F)
pathway=as.matrix(pathway)
pathway=pathway[-c(20:24,29,32,34,63,64,78,170)]
pw1=paste("extdata",pathway,sep="/")
pw2=paste(pw1,"xml",sep=".")
I=as.character(1:length(pw2))
len=0
for (i in 1:length(pw2)){
  #sfile <- system.file(pw2[i],package="KEGGgraph") # OLD
  sfile <- pw2[i]
  gR <- parseKGML2Graph(sfile,expandGenes=TRUE)
  len[i]=length(nodes(gR))
}

pathway1=matrix(c(pathway,len),ncol=2)
subpathway1=subset(pathway1,as.numeric(pathway1[,2])<1000)
subpathway2=subset(subpathway1,as.numeric(subpathway1[,2])>10)
maxlen=max(as.numeric(subpathway2[,2]))
subpathway=subpathway2[,-2]
S1=matrix(,nrow=length(subpathway),ncol=maxlen,dimname=list(subpathway))
subpw1=paste("extdata",subpathway,sep="/")
subpw2=paste(subpw1,"xml",sep=".")
for (i in 1:length(subpathway)){
  # sfile <- system.file(subpw2[i],package="KEGGgraph") # OLD
  sfile <- subpw2[i]
  gR <- parseKGML2Graph(sfile,expandGenes=TRUE)
  ge=nodes(gR)
  S1[i,1:length(ge)]=t(ge)
}

len1=len[which(len<1000)]
len2=len1[which(len1>10)]

###############################################;

dat1=read.csv(file="dat.csv", header=F)
cn=paste("hsa:",dat1[1,],sep="")
dat1=as.matrix(dat1)
dat2=dat1[-1,]
colnames(dat2)=cn

library(impute)

dat3=impute.knn(t(dat2))
dat4=t(dat3$data)
Y=c(4.6, 3.9, 3.9, 6.2, 4.6, 4.6, 4.8, 13.6, 9.8, 4.6, 2.6, 7.1, 6.8, 5.4, 0, 9.3, 0.09, 3, 7, 9.8, 4.7, 3.5, 0, 5.4, 44.6, 22, 5.2, 7, 7.8, 4.6, 8.5, 4.7, 17.8, 0.3, 24.1, 12.4, 4.2, 4.3, 13.6, 5.4, 6, 6.8, 9.5, 2.6, 5.6, 3.8, 0, 26.8, 6.4, 14.9, 9.7, 3.9, 6.9, 5.9, 11.4, 7.8, 5.8, 14.7, 8.5, 25, 33.4, 23.9, 16.5, 15.5, 57.5, 9.6, 0, 57, 48.4, 38, 43.8, 15, 12.1, 26.2, 12.1, 5, 5.7, 7.3, 5, 5, 1.7)
which(Y==0)
dat5=dat4[-c(15,23,47,67,82:101),]
Y=Y[which(Y!=0)]



dat5[,cn[match(S1[1,1:len2[1]],cn,nomatch=0)]]


X1=dat5[,cn[match(S1[1,1:len2[1]],cn,nomatch=0)]]
X2=dat5[,cn[match(S1[2,1:len2[2]],cn,nomatch=0)]]

MATS<-list(X1,X2)

for(i in 3:219){
  X3=dat5[,cn[match(S1[i,1:len2[i]],cn,nomatch=0)]]
  MATS[[i]]<-X3
}

p.data=dim(MATS[[1]])[2]
a_init = 0.1;
b_init = 0.1;
r_init = 0.1;
ahyper=0.001
bhyper=0.001
rhyper=0.001
Beta_init = rnorm(p.data,0,0);
sigma2_init = var(Y);
Lambda_init = matrix(0, p.data, p.data); 
diag(Lambda_init) = 1;
Nburn = 10000;
Niter = 20000;
Nchain = 1;
Nthin=1;

MCMCResult = BVSGR(MATS[[1]], Y, a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin);
betahat1=MCMCResult[[1]]$Betahat
betahat1
myheatmap(MATS[[1]],MCMCResult[[1]]$Lambdahat)
Lambda1=MCMCResult[[1]]$Lambdahat
simga2=mean(MCMCResult[[1]]$sigma2hist)

p.data=dim(MATS[[2]])[2]
Beta_init = rnorm(p.data,0,0);
sigma2_init = var(Y);
Lambda_init = matrix(0, p.data, p.data); 
diag(Lambda_init) = 1;
MCMCResult = BVSGR(MATS[[2]], Y, a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin);
betahat2=MCMCResult[[1]]$Betahat
betahat2
myheatmap(MATS[[2]],MCMCResult[[1]]$Lambdahat)
Lambda2=MCMCResult[[1]]$Lambdahat
simga2[2]=mean(MCMCResult[[1]]$sigma2hist)


BETA=list(betahat1,betahat2)
LAMBDA=list(Lambda1,Lambda2)
for(i in 3:219){
  p.data=dim(MATS[[i]])[2]
  Beta_init = rnorm(p.data,0,0);
  sigma2_init = var(Y);
  Lambda_init = matrix(0, p.data, p.data); 
  diag(Lambda_init) = 1;
  MCMCResult = BVSGR(MATS[[i]], Y, a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin);
  BETA[[i]]=betahat=MCMCResult[[1]]$Betahat
  
  LAMBDA[[i]]=MCMCResult[[1]]$Lambdahat
  simga2[i]=mean(MCMCResult[[1]]$sigma2hist)
}

p=193
BETA[[p]]
pdf("mygraph.pdf")
for (p in 1:219){
  myheatmap(MATS[[p]],LAMBDA[[p]])
}
dev.off()


ipath=c(81, 82, 99, 103, 108, 111, 134, 136, 145, 193, 195, 198, 199, 200, 205, 208)
thhld=seq(.1,.2,.01)

for(i in 1:11){
  label_vector <-c(seq(1:dim( LAMBDA[[ipath[1]]])[1] ))
  mymat      <-solve(t(MATS[[ipath[1]]] )%*%MATS[[ipath[1]]] + LAMBDA[[ipath[1]]])
  mysds      <-sqrt(diag(mymat))
  mycorrmat  <-mymat/outer(mysds,mysds,FUN="*")
  
  source("myImagePlot.r")
  
  myImagePlot(abs(mycorrmat))
  
  gplot(abs(mycorrmat),edge.lwd= abs(mycorrmat*20),gmode="graph",pad=0.3,loop.cex=5,thresh=thhld[i],displaylabels=TRUE,label=label_vector,label.bg ="gray90",edge.col="blue",xlab="Dependence Structure among Covariates")
  
}

p=16
subpathway2[ipath[p],]
myheatmap(MATS[[ipath[p]]],LAMBDA[[ipath[p]]])
#BETA[[ipath[p]]]
colnames(MATS[[ipath[p]]])[which(abs(BETA[[ipath[p]]])>2)]

p.data=dim(MATS[[ipath[p]]])[2]
Beta_init = rnorm(p.data,0,0);
sigma2_init = var(Y);
Lambda_init = matrix(0, p.data, p.data); 
diag(Lambda_init) = 1;
MCMCResult = BVSGR(MATS[[ipath[p]]], Y, a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin);
betahat1=MCMCResult[[1]]$Betahat
betahat1
myheatmap(MATS[[ipath[p]]], MCMCResult[[1]]$Lambdahat)
Lambda1=MCMCResult[[1]]$Lambdahat
#sum(MCMCResult[[1]]$Betahist[,54]>0)/10000



######################

p=2
label_vector <-c(seq(1:dim( LAMBDA[[ipath[p]]])[1] ))
mymat      <-solve(t(MATS[[ipath[p]]] )%*%MATS[[ipath[p]]] + LAMBDA[[ipath[p]]])
mysds      <-sqrt(diag(mymat))
mycorrmat  <-mymat/outer(mysds,mysds,FUN="*")

source("myImagePlot.r")

myImagePlot(abs(mycorrmat))

gplot(abs(mycorrmat),edge.lwd= abs(mycorrmat*20),gmode="graph",pad=0.3,loop.cex=5,thresh=.15,displaylabels=TRUE,label=label_vector,label.bg ="gray90",edge.col="blue",xlab="Dependence Structure among Covariates")

mat=mycorrmat
mat1=lower.tri(mat, diag=T)
mat[mat1]<-NA
loc1=which(abs(mat)>.15, arr.ind=T)
genename=unique(c(rownames(mat)[loc1[,1]],colnames(mat)[loc1[,2]]))
genename


######################################################

set.seed(1)
training=sample(77,57)
p.data=dim(MATS[[ipath[1]]])[2]
a_init = 0.1;
b_init = 0.1;
r_init = 0.1;
ahyper=0.001
bhyper=0.001
rhyper=0.001
Beta_init = rnorm(p.data,0,0);
sigma2_init = var(Y);
Lambda_init = matrix(0, p.data, p.data); 
diag(Lambda_init) = 1;
Nburn = 10000;
Niter = 20000;

MCMCResult = BVSGR(MATS[[ipath[1]]][training,], Y[training], a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin);
betahat1=MCMCResult[[1]]$Betahat
betahat1
MSE=sum(MATS[[ipath[1]]][-training,]%*%betahat1-Y[-training])^2/20

MSE=matrix(,nrow=50,ncol=16)
set.seed(1)
for (t in 1:50){
  training=sample(77,57)
  for(i in 1:16){
    p.data=dim(MATS[[ipath[i]]])[2]
    a_init = 0.1;
    b_init = 0.1;
    r_init = 0.1;
    ahyper=0.001
    bhyper=0.001
    rhyper=0.001
    Beta_init = rnorm(p.data,0,0);
    sigma2_init = var(Y);
    Lambda_init = matrix(0, p.data, p.data); 
    diag(Lambda_init) = 1;
    Nburn = 10000;
    Niter = 20000;
    X=scale(MATS[[ipath[i]]])
    Y1=scale(Y,scale=F)
    MCMCResult = BVSGR(X[training,], Y1[training], a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin);
    betahat1=MCMCResult[[1]]$Betahat
    betahat1
    MSE[t,i]=sum(X[-training,]%*%betahat1-Y1[-training])^2/20
  }
}

set.seed(1)
training=sample(77,57)
for(i in 1:16){
  p.data=dim(MATS[[ipath[i]]])[2]
  a_init = 0.1;
  b_init = 0.1;
  r_init = 0.1;
  ahyper=0.001
  bhyper=0.001
  rhyper=0.001
  Beta_init = rnorm(p.data,0,0);
  sigma2_init = var(Y);
  Lambda_init = matrix(0, p.data, p.data); 
  diag(Lambda_init) = 1;
  Nburn = 10000;
  Niter = 20000;
  
  MCMCResult = BVSGR(MATS[[ipath[i]]][training,], Y[training], a_init, b_init, r_init, Beta_init, Lambda_init, sigma2_init,ahyper,bhyper,rhyper, Nburn, Niter,Nchain,Nthin);
  betahat1=MCMCResult[[1]]$Betahat
  betahat1
  MSE[i]=sum(MATS[[ipath[i]]][-training,]%*%betahat1-Y[-training])^2/20
}

#####################lasso
library(lars)
lasso=lars(MATS[[ipath[i]]][training,],Y[training],intercept=F)
X=scale(MATS[[ipath[i]]][training,])
Y1=scale(Y[training],scale=F)
lasso=lars(X,Y1,intercept=F,normalize=F)

MSEls=matrix(,nrow=50,ncol=16)
set.seed(1)
for (t in 1:50){
  training=sample(77,57)
  for(i in 1:16){

    X=scale(MATS[[ipath[i]]])
    Y1=scale(Y,scale=F)
    lasso = lars(X[training,], Y1[training],intercept=F,normalize=F);
    betahat1=lasso$beta[dim(lasso$beta)[1],]
    MSEls[t,i]=sum(X[-training,]%*%betahat1-Y1[-training])^2/20
  }
}



###########Blasso
library(monomvn)


MSE=matrix(,nrow=50,ncol=16)
set.seed(1)
for (t in 1:50){
  training=sample(77,57)
  for(i in 1:16){
    blassoresult=blasso(MATS[[ipath[i]]][training,],Y[training], T=5000)
    
    betahat1=colMeans(blassoresult$beta[2501:5000,])
    mu1=mean(blassoresult$mu[2501:5000])
    MSE[t,i]=sum(MATS[[ipath[i]]][-training,]%*%betahat1+mu1-Y[-training])^2/20
  }
}
