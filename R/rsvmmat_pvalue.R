#' Calculate prospective and retrospective P-values 
#'
#' This function tests a SNPs for a given SNP set for a given rsvmmat estimated null model.
#'
#' @param rsvmmat.est The output of function "rsvmmat_est()"
#' @param G The genotype matrix, an m*q matrix where m is the number of subjects and q is the total number genetic variants. 
#' @param impute.method choose the iputation method when there is missing genotype. Optional options are: 'random', 'fixed' or 'bestguess'.
#' @param GRM takes m-by-m genetic correlation matrix or kinship matrix.
#' 
#' @return This function returns a vector of the pvalue of SVMMAT-B, SVMMAT-S,SVMMAT-A, RSVMMAT-B, RSVMMAT-S, and RSVMMAT-A.
#' 
#' @export

rsvmmat_test <-function(rsvmmat.est, G, impute.method='fixed', GRM = NULL)
{
  res<-rsvmmat.est$Y.res; phi=rsvmmat.est$phi; X<-rsvmmat.est$X;N<-nrow(X)
  m<-rsvmmat.est$m;time<-rsvmmat.est$time;mu<-rsvmmat.est$mu;tau<-rsvmmat.est$tau;cluster.id<-rsvmmat.est$cluster.id
  snp.names<-colnames(G); family = rsvmmat.est$family;Rr=rsvmmat.est$Rr;Td=rsvmmat.est$Td;Hd=rsvmmat.est$Hd
  ntime=length(unique(time[,2]))
  
  if(!is.null(GRM)){
    dGRM=sqrt(diag(GRM))
    for(i in 1:dim(GRM)[1]) {
      GRM[,i]=GRM[,i]/dGRM[i]
      GRM[i,]=GRM[i,]/dGRM[i]
    }
    
  }
  
  
  G[G==9]<-NA
  N_MISS<-sum(is.na(G))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G<-Impute(G,'fixed')
  }
  
  var.G<-t(t(G)-colMeans(G));
  center.G = var.G
  SNP.list<-colnames(G)
  MAF<-colMeans(as.matrix(G[,colMeans(var.G^2)!=0]), na.rm=T)/2;
  G<-as.matrix(G[,colMeans(var.G^2)!=0]); center.G = center.G[,colMeans(var.G^2)!=0]
  
  SNP.name<-SNP.list[colMeans(center.G^2)!=0]
  
  weights = rep(0,ncol(G));
  index.common = MAF>=0.05
  index.rare = MAF<0.05
  scaler = stats::dbeta(0.05,1,25)/stats::dbeta(0.05,0.5,0.5)
  weights[index.common] = scaler*stats::dbeta(MAF[index.common],0.5,0.5)
  weights[index.rare] = stats::dbeta(MAF[index.rare],1,25)
  
  #weights<-dbeta(MAF,1,25); 
  #WG=G%*%weights;WG<-WG[match(time[,1],cluster.id)]
  WG<-t(t(G)*weights)
  var_g <-stats::var(WG%*%rep(1,dim(G)[2]))[1,1]
  Var_WG=stats::cov(WG)
  
  WG<-WG[match(time[,1],cluster.id),]
  P1<-rsvmmat.est$P1; P2<-rsvmmat.est$P2;
  # Z<-G-P2%*%(P1%*%G);
  
  ResN=c(res)*rsvmmat.est$inciN;  ResNTd=ResN%*%Td
  BTResNTd=matrix(0,nrow=m,ncol=Hd)
  
  IDtable=table(time[,1])
  for(ss in 1:m){
    if(IDtable[ss]==1){
      BTResNTd[ss,]= ResNTd[which(time[,1]==ss),]
    }else{ BTResNTd[ss,]= colSums(ResNTd[which(time[,1]==ss),]) 
    }
  }
  
  halfR=Re(expm::sqrtm(Rr)); ResNR=ResN%*%halfR; 
  BTResNR=matrix(rep(0,m*ntime),ncol=ntime)
  
  for(ss in 1:m){
    if(IDtable[ss]==1){
      BTResNR[ss,]= ResNR[which(time[,1]==ss),]
    }else{ BTResNR[ss,]= colSums(ResNR[which(time[,1]==ss),]) 
    }
  }
  
  BTResNR2colS=colSums(BTResNR^2)
  
  
  if(is.null(GRM)){Var.retro=t(BTResNTd)%*%BTResNTd
  V.Rretro=t(BTResNR)%*%BTResNR
  }else{       Var.retro=t(BTResNTd)%*%GRM%*%BTResNTd
  V.Rretro=t(BTResNR)%*%GRM%*%BTResNR}
  
  type1result=NULL
  GN=c(WG%*%rep(1,dim(G)[2]))*rsvmmat.est$inciN
  GNTd=GN%*%Td; colnames(GNTd)=1:Hd;tranG=GN%*%halfR;GNR=GN%*%Rr;
  ZGNTd=GNTd- P2%*%(P1%*%GNTd); ZGNR=tranG-P2%*%(P1%*%tranG)
  
  V.pro1<-t(ZGNTd)%*%rsvmmat.est$V.inv%*%ZGNTd;
  V.RrB=t(ZGNR)%*%rsvmmat.est$V.inv%*%ZGNR;   
  tran_res<-res
  
  score=(t(tran_res)%*%GNTd%*%MASS::ginv(V.pro1)%*%t(GNTd)%*%tran_res)[1,1]/(phi^2)
  score.retro=(t(tran_res)%*%GNTd%*%MASS::ginv(Var.retro*var_g)%*%t(GNTd)%*%tran_res)[1,1]
  pHotellingsT=stats::pchisq(score,df=Hd,lower.tail=F)
  pRHotellingsT=stats::pchisq(score.retro,df=Hd,lower.tail=F)      
  
  HalfTBurden=t(tranG)%*%tran_res
  TBurden=sum(HalfTBurden^2)
  lambdaS=eigen(V.RrB,symmetric = TRUE,only.values = TRUE)$values
  
  if(sum(lambdaS)==0) { pSmoothsSKAT2=1
  }else{ pTBurden=generalchisq(lambdaS,TBurden/(phi^2))}
  
  lambdaS.retro=eigen(V.Rretro*var_g,symmetric = TRUE,only.values = TRUE)$values
  if(sum( lambdaS.retro)==0) {pSmoothsSKAT.retro2=1
  }else{
    pTrBurden=generalchisq(lambdaS.retro,TBurden)}
  
  Tp=(tan((0.5-pHotellingsT)*pi)+tan((0.5-pTBurden)*pi))/2
  pSVMMATB=0.5-atan(Tp)/pi
  Tr=(tan((0.5-pRHotellingsT)*pi)+tan((0.5-pTrBurden)*pi))/2
  pSRVMMATB=0.5-atan(Tr)/pi
  
  
  TSKAT=0;ALLZGNR=NULL
  for(k in 1:dim(G)[2]){
    GN=WG[,k]*rsvmmat.est$inciN;tranG=GN%*%halfR;ZGNR=tranG-P2%*%(P1%*%tranG)
    # V.Rr=t(ZGNR)%*%V.inv%*%ZGNR;   
    HalfTSKAT=t(tranG)%*%tran_res
    TSKAT=TSKAT+sum(HalfTSKAT^2)
    ALLZGNR=cbind(ALLZGNR,ZGNR)
  }
  ALLZGNRV=t(ALLZGNR)%*%rsvmmat.est$V.inv
  V.Rr=ALLZGNRV%*%ALLZGNR
  lambdaS=eigen(V.Rr,symmetric = TRUE,only.values = TRUE)$values
  
  if(sum(lambdaS)==0) { pSmoothsSKAT2=1
  }else{ pTSKAT=generalchisq(lambdaS,TSKAT/(phi^2))}
  
  lambdaA=eigen(Var_WG,symmetric = TRUE,only.values = TRUE)$values
  lambdaB=eigen(V.Rretro,symmetric = TRUE,only.values = TRUE)$values
  lambdaS.retro=kronecker(lambdaA,lambdaB)
  if(sum( lambdaS.retro)==0) {pSmoothsSKAT.retro2=1
  }else{
    pTrSKAT=generalchisq(lambdaS.retro,TSKAT)}
  
  Tp=(tan((0.5-pHotellingsT)*pi)+tan((0.5-pTSKAT)*pi))/2
  pSVMMATS=0.5-atan(Tp)/pi
  Tr=(tan((0.5-pRHotellingsT)*pi)+tan((0.5-pTrSKAT)*pi))/2
  pSRVMMATS=0.5-atan(Tr)/pi
  
  Tp=(tan((0.5-pSVMMATB)*pi)+tan((0.5-pSVMMATS)*pi))/2
  pSVMMATA=0.5-atan(Tp)/pi
  Tr=(tan((0.5-pSRVMMATB)*pi)+tan((0.5-pSRVMMATS)*pi))/2
  pSRVMMATA=0.5-atan(Tr)/pi
  
  
  type1result=cbind(pSVMMATB,pSVMMATS,pSVMMATA,pSRVMMATB,pSRVMMATS,pSRVMMATA)
  colnames(type1result)=c("SVMMATB","SVMMATS","SVMMATA","RSVMMATB","RSVMMATS","RSVMMATA")
  return(type1result)
}




generalchisq=function(lambda,Q){
  muq=sum(lambda)
  sigmaq=sqrt(2*sum(lambda^2))
  s1=sum(lambda^3)/sum(lambda^2)^1.5
  s2=sum(lambda^4)/sum(lambda^2)^2
  if((s1^2)>s2){
    a=1/(s1-sqrt(s1^2-s2))
    delta=s1*a^3-a^2
    l=a^2-2*delta
  }else{
    delta=0
    l=sum(lambda^2)^3/sum(lambda^3)^2
  }
  stats::pchisq((sum(Q)-muq)/sigmaq*sqrt(2*(l+2*delta))+l+delta,df=l, ncp=delta,lower.tail=FALSE)
}





Impute<-function(Z, impute.method){
 if(is.vector(Z)){
      if(impute.method =="random"){
            IDX<-which(is.na(Z))
            if(length(IDX) > 0){
                maf1<-mean(Z[-IDX])/2
                Z[IDX]<-stats::rbinom(length(IDX),2,maf1)
            }
      } else if(impute.method =="fixed"){
              IDX<-which(is.na(Z))
              if(length(IDX) > 0){
                   maf1<-mean(Z[-IDX])/2
                   Z[IDX]<-2 * maf1
              }
       } else if(impute.method =="bestguess") {
               IDX<-which(is.na(Z))
               if(length(IDX) > 0){
                    maf1<-mean(Z[-IDX])/2
                   Z[IDX]<-round(2 * maf1)
               }
       } else {
             stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
       }
  } else{
      p<-dim(Z)[2]
      if(impute.method =="random"){
         for(i in 1:p){
            IDX<-which(is.na(Z[,i]))
            if(length(IDX) > 0){
                maf1<-mean(Z[-IDX,i])/2
                Z[IDX,i]<-stats::rbinom(length(IDX),2,maf1)
            }
         }
      } else if(impute.method =="fixed"){
          for(i in 1:p){
              IDX<-which(is.na(Z[,i]))
              if(length(IDX) > 0){
                   maf1<-mean(Z[-IDX,i])/2
                   Z[IDX,i]<-2 * maf1
              }
          }
       } else if(impute.method =="bestguess") {
           for(i in 1:p){
               IDX<-which(is.na(Z[,i]))
               if(length(IDX) > 0){
                    maf1<-mean(Z[-IDX,i])/2
                   Z[IDX,i]<-round(2 * maf1)
               }
          }
       } else {
             stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
       }
  }
  return(as.matrix(Z))
}

