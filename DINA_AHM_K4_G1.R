library("Rcpp", lib.loc="~/R_libs")
library("RcppArmadillo", lib.loc="~/R_libs")
library("gtools", lib.loc="~/R_libs")


Rcpp::sourceCpp('./DINA_AHM_v10.cpp')
args <- commandArgs(trailingOnly = TRUE)
rep = as.integer(args)


N=500
K=4
J=20

vec=c(1,1,0,0,
        1,0,1,0,
        1,0,0,1,
        0,1,1,0,
        0,1,0,1,
        0,0,1,1,
        1,1,1,0,
        1,1,0,1,
        1,0,1,1,
        0,1,1,1,
        1,1,1,1,
        1,1,1,1)
Q1=matrix(vec,ncol=4,byrow=TRUE)
Q=rbind(diag(4),diag(4),Q1)

G=matrix(0,ncol=K,nrow=K)
G1=G;G2=G;G3=G;G4=G;
G1[1,2]=1;G1[2,3]=1;G1[3,4]=1
G2[1,2]=1;G2[2,4]=1;G2[1,3]=1;G2[3,4]=1
G3[1,2]=1;G3[1,3]=1;G3[1,4]=1
G4[1,2]=1;G4[1,3]=1;G4[3,4]=1;

ss=rep(0.2,J)
gs=rep(0.2,J)


#
G0=G1

R0=Reachability(G0,K)
C0=ConnectMat(R0,K)

rho=0
Z = matrix(rnorm(N*K),N,K)
Sig = matrix(rho,K,K)
diag(Sig) = 1
theta = Z%*%chol(Sig)
thvals = 0#matrix(qnorm((1:K)/(K+1)),N,K,byrow=T)
Alphas = 1*(theta>thvals)
vv=bijectionvector(K)
PIndex=Collapsed(R0,K)
CLASS01=Alphas%*%vv
CLASS0=PIndex$NewInd[CLASS01+1]
PI0=table(CLASS0)/N
PI0
ETA=ETAmat(K,J,Q)

Y=sim_Y_dina(N,J,CLASS0,ETA,gs,ss)


p1=0.3;p2=0.3;
ptm <- proc.time()
Out=dina_AHM(Y,Q,4,p1,p2,15000,20000)
proc.time() - ptm

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

table(Out$Trace[1,])
table(Out$Trace[2,])
table(Out$Trace[2,])/table(Out$Trace[1,])

gs_est=apply(Out$GS,1,mean);gs_est
ss_est=apply(Out$SS,1,mean);ss_est
G_est=apply(Out$Graphs,c(1,2),mean);G_est
R_est=apply(Out$RS,c(1,2),mean);R_est
PIs_est=apply(Out$PIs,1,mean);PIs_est

Class_est<-apply(Out$CLASS,1,getmode)
PAR=sum(Class_est==CLASS0)/N;PAR
Alphas_est<-sapply(Class_est,inv_bijectionvector, K = K)
Alphas0<-sapply(CLASS0,inv_bijectionvector, K = K)
AAR=apply(Alphas_est==Alphas0,1,mean);AAR
SS_RMSE=sqrt(mean((ss_est-ss)^2))
GS_RMSE=sqrt(mean((gs_est-gs)^2))


sink(paste("./G1rho0","N=",N, "K=",K,
           "rep",rep, ".txt",sep=''))

cat("Estimated G")
print(G_est)
cat("\n")
cat("Estimated R")
print(R_est)
cat("\n")
cat("PAR\n")
print(PAR)
cat("\n")
cat("AAR\n")
print(AAR)
cat("\n")
cat("DINA parameter ss\n")
print(ss_est)
cat("\n")
cat("DINA parameter gs\n")
print(gs_est)
cat("\n")
cat("PI True\n")
print(PI0)
cat("\n")
cat("PI\n")
print(PIs_est)
cat("\n")
cat("ss_RMSE\n")
print(SS_RMSE)
cat("\n")
cat("GS_RMSE\n")
print(GS_RMSE)
cat("\n")
sink()
###################################
