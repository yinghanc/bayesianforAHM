#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
double MAX=1e300;
double MIN=1e-300;
// [[Rcpp::export]]
arma::vec bijectionvector(unsigned int K) {
  arma::vec vv(K);
  for(unsigned int k=0;k<K;k++){
    vv(k) = pow(2,K-k-1);
  }
  return vv;
}

// [[Rcpp::export]]
arma::vec inv_bijectionvector(unsigned int K,double CL){
  arma::vec alpha(K);
  for(unsigned int k=0;k<K;k++){
    double twopow = pow(2,K-k-1);
    alpha(k) = (twopow<=CL);
    CL = CL - twopow*alpha(k);
  }
  return alpha;
}

// [[Rcpp::export]]
arma::mat inv_bijectionmat(unsigned int K,const arma::vec& CL) {
  unsigned int Col=CL.n_elem;
  arma::mat alpha(K,Col);
  for(unsigned int i=0;i<Col;i++){
    arma::colvec alphai(K);double cl=CL(i);
    for(unsigned int k=0;k<K;k++){
      double twopow = pow(2,K-k-1);
      alphai(k) = (twopow<=cl);
      cl = cl - twopow*alphai(k);
    }
    alpha.col(i)=alphai;
  }
  return alpha;
}

// [[Rcpp::export]]
arma::mat ETAmat(unsigned int K,unsigned int J,const arma::mat& Q) {
  double nClass = pow(2,K);
  arma::mat ETA(J,nClass);
  for(unsigned int cc=0;cc<nClass;cc++){
    arma::vec alpha_c = inv_bijectionvector(K,cc);
    
    for(unsigned int j=0;j<J;j++){
      arma::rowvec qj = Q.row(j);
      double compare = arma::conv_to<double>::from(qj*alpha_c - qj*qj.t() ); 
      ETA(j,cc) = (compare>=0);
    }
  }
  return ETA;
}

// [[Rcpp::export]]
arma::vec abcounts(unsigned int N,const arma::vec& Yj,
                   const arma::vec& CLASS,const arma::vec& ETAtnokimes1ma){
  arma::vec ab = arma::zeros<arma::vec>(2);
  for(unsigned int i=0;i<N;i++){
    ab(Yj(i)) += ETAtnokimes1ma(CLASS(i));
  }
  return ab;
}

// [[Rcpp::export]]
arma::mat ClassbyQmat(unsigned int K){
  double nClass = pow(2,K);
  arma::mat ETAbyQ(nClass,nClass-1);
  for(unsigned int cc=0;cc<nClass;cc++){
    arma::vec alpha_c = inv_bijectionvector(K,cc);
    for(unsigned int r=0;r<nClass-1;r++){
      arma::vec qj = inv_bijectionvector(K,r+1);
      double compare = arma::conv_to<double>::from(qj.t()*alpha_c - qj.t()*qj ); 
      ETAbyQ(cc,r) = (compare>=0);
    }
  }
  return ETAbyQ;
}

// [[Rcpp::export]]
double rmultinomial(const arma::vec& ps){
  unsigned int C = ps.n_elem;
  double u = R::runif(0,1);
  arma::vec cps = cumsum(ps);
  arma::vec Ips = arma::zeros<arma::vec>(C);
  Ips.elem(arma::find(cps < u) ).fill(1.0);
  
  return sum(Ips);
}

// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec& deltas){
  unsigned int C = deltas.n_elem;
  arma::vec Xgamma(C);
  
  //generating gamma(deltac,1)
  for(unsigned int c=0;c<C;c++){
    Xgamma(c) = R::rgamma(deltas(c),1.0);
  }
  return Xgamma/sum(Xgamma);
}

// [[Rcpp::export]]
double pYi(const arma::vec& ETA_i,const arma::vec& Y_i,const arma::vec& ss,
            const arma::vec& gs){
  arma::vec one_m_ss = 1. - ss;
  arma::vec one_m_gs = 1. - gs;
  arma::vec one_m_ETA_i = 1. - ETA_i;
  arma::vec one_m_Y_i = 1. - Y_i;
  
  arma::vec ps = Y_i%(one_m_ss%ETA_i + gs%one_m_ETA_i) + one_m_Y_i%(ss%ETA_i + one_m_gs%one_m_ETA_i);
  
  return arma::prod(ps);
}

// [[Rcpp::export]]
double pYiratio(const arma::vec& ETA1,const arma::vec& ETA2,const arma::vec& Y_i,const arma::vec& ss,
           const arma::vec& gs){
  arma::uvec ind = find(ETA1!=ETA2);
  arma::vec one_m_ss = 1. - ss(ind);
  arma::vec one_m_gs = 1. - gs(ind);
  arma::vec one_m_ETA1 = 1. - ETA1(ind);
  arma::vec one_m_ETA2 = 1. - ETA2(ind);
  arma::vec one_m_Y_i = 1. - Y_i(ind);
  
  double ratio=1.0;
  if(ind.n_elem>0){
    arma::vec ps1 = Y_i(ind)%(one_m_ss%ETA1(ind) + gs(ind)%one_m_ETA1) + one_m_Y_i%(ss(ind)%ETA1(ind) + one_m_gs%one_m_ETA1);
    arma::vec ps2 = Y_i(ind)%(one_m_ss%ETA2(ind) + gs(ind)%one_m_ETA2) + one_m_Y_i%(ss(ind)%ETA2(ind) + one_m_gs%one_m_ETA2);
    ratio=arma::prod(ps1)/arma::prod(ps2);
  }
  return ratio;
}

// [[Rcpp::export]]
arma::mat random_Q(unsigned int J,unsigned int K){
  
  //Generate identity matrices
  arma::vec one_K = arma::ones<arma::vec>(K);
  arma::mat I_K = arma::diagmat(one_K);
  //arma::mat Two_I_K = arma::join_cols(I_K,I_K);
  
  //generate Q1
  unsigned int JmK = J-K;
  unsigned int J1max = K;
  if(K>JmK){
    J1max = JmK;
  }
  unsigned int J1 = arma::conv_to< unsigned int >::from(arma::randi<arma::vec>(1,arma::distr_param(1,J1max) ) );
  arma::mat U1 = arma::randu<arma::mat>(J1,K);
  arma::mat Q1 = arma::zeros<arma::mat>(J1,K);
  
  //fix elements so columns are nonzero
  arma::vec row_ks = arma::randi<arma::vec>(K,arma::distr_param(0,J1-1) );
  for(unsigned int k=0;k<K;k++){
    Q1(row_ks(k),k) = 1;
  }
  
  Q1.elem(arma::find(Q1 > .5) ).fill(1.0);
  
  arma::mat Q = arma::join_cols(I_K,Q1);
  
  //Generating the remaining elements of Q in Q2 
  unsigned int JmKmJ1 = JmK - J1;
  arma::mat Q2 = arma::zeros<arma::mat>(JmKmJ1,K);
  if(JmKmJ1>0){
    arma::mat U2 = arma::randu<arma::mat>(JmKmJ1,K);
    Q2.elem(arma::find(U2 > .5) ).fill(1.0);
    Q = arma::join_cols(Q,Q2);
  }
  
  //Q
  arma::uvec P = arma::uvec(J);
  for(unsigned int j=0;j<J;j++){
    P(j)=j;
  }
  P = arma::shuffle(P);
  return Q.rows(P);
}

// [[Rcpp::export]]
arma::mat sim_Y_dina(unsigned int N,unsigned int J,const arma::vec& CLASS,
                   const arma::mat& ETA,const arma::vec& gs,const arma::vec& ss){
  arma::mat Y(N,J);
  for(unsigned int i=0;i<N;i++){
    double class_i = CLASS(i);
    arma::vec ETA_i = ETA.col(class_i);
    for(unsigned int j=0;j<J;j++){
      double u = R::runif(0,1);
      Y(i,j) = 1.*(gs(j)*(1.-ETA_i(j)) + (1.-ss(j))*ETA_i(j) > u);
    }
  }
  return Y;
}


// [[Rcpp::export]]
double dina_m2LL(unsigned int N, unsigned int J, unsigned int K, const arma::mat &Y, const arma::vec &pis,
                 const arma::mat Q, const arma::vec& ss, const arma::vec& gs)
{ double m2ll = 0.;
  unsigned int nClass=pow(2, K);
  arma::mat ETA = ETAmat(K,J,Q);
  arma::vec one_m_ss = 1. - ss;
  arma::vec one_m_gs = 1. - gs;
  for (unsigned int i = 0; i < N; i++) {
    arma::vec Yi =(Y.row(i)).t();
    double pYi = 0.;
    for (unsigned int cc = 0; cc < nClass; cc++) {
      double py_c = 1.;
      arma::vec ETA_c = ETA.col(cc);
      arma::vec one_m_ETA_c = 1. - ETA_c;
      arma::vec one_m_Yi = 1. - Yi;
      arma::vec ps = Yi%(one_m_ss%ETA_c + gs%one_m_ETA_c) + one_m_Yi%(ss%ETA_c + one_m_gs%one_m_ETA_c);
      py_c=arma::prod(ps);
      pYi += py_c * pis(cc);
    }
    m2ll += log(pYi);
  }
  return -2. * m2ll;
}

// [[Rcpp::export]]
arma::mat Boolean(arma::mat A){
  arma::mat R;
  R=arma::conv_to<arma::mat>::from(A>0);
  return R;
}

// [[Rcpp::export]]
arma::vec Booleanvec(arma::vec A){
  arma::vec R;
  R=arma::conv_to<arma::vec>::from(A>0);
  return R;
}

// [[Rcpp::export]]
arma::mat Reachability(const arma::mat& StrucMat,unsigned int K){
  arma::mat R=StrucMat;
    R.diag().zeros();
    arma::mat Identity=arma::eye(K,K);
    arma::mat Rnext=R+Identity;
    while(accu(R==Rnext)<(K*K)){
      R=Rnext;
      Rnext=Boolean(R*(R+Identity));
    }
    R.diag().ones();
  return R;
}


// [[Rcpp::export]]
arma::mat ConnectMat(const arma::mat& R,unsigned int K){
  arma::mat Connect(K,K);Connect.eye();
  for(unsigned int i=0;i<(K-1);i++){
    for(unsigned int j=(i+1);j<K;j++){
      if((R(i,j)==1)||(R(j,i)==1)){
        Connect(i,j)=1;Connect(j,i)=1;
      }
    }
  }
  return Connect;
}


// [[Rcpp::export]]
int check_tri(const arma::mat& StrucMat,unsigned int K){//check if triangle free
  int flag=0;
  arma::mat temp(K,K);temp.zeros();
  for(unsigned int i=0;i<K;i++){
    for(unsigned int j=0;j<K;j++){
      if(StrucMat(i,j)>0){
        temp(i,j)=StrucMat(i,j);
        temp(j,i)=StrucMat(i,j);
      }
    }
  }
  temp.diag().zeros();
  double tr;
  tr=trace(temp*temp*temp);
  if(tr==0){flag=1;}else{flag=0;}
  return flag;
}

// [[Rcpp::export]]
int check_AHM(const arma::mat& R,unsigned int K){//Hierarchy structure: sibling cannot be parents 
  int flag=1;
  arma::mat temp=R;temp.diag().zeros();
  for(unsigned int i=0;i<K;i++){
    arma::rowvec Rrowi=temp.row(i);
    arma::uvec index=find(Rrowi>0);
    if(index.n_elem==0){
      continue;
    }else{
      arma::vec Rcoli=temp.col(i);
      double check=accu(Rcoli(index));
      if(check>0){
        flag=0;return flag;
      }else{
        continue;
      }
    }
  }
  return flag;
}


 // [[Rcpp::export]]
 int check_G(const arma::mat& StrucMat,const arma::mat& R,unsigned int K){
 int flag=1;
 //int check0=check_1D(StrucMat,K);
 //if(check0<1){flag=0;return flag;}
// int check1=check_tri(StrucMat,K);
// if(check1<1){flag=0;return flag;}
 int check2=check_AHM(R,K);
 if(check2<1){flag=0;return flag;}
// if((check1+check2)==2){flag=1;}else{flag=0;}
 return flag;
 }


// [[Rcpp::export]]
arma::uvec check_no_connect(const arma::mat& C,unsigned int K){
  arma::uvec q=find(C==0);
  return q;
}

// [[Rcpp::export]]
arma::uvec check_link(const arma::mat& G,unsigned int K){
  arma::mat ModG=G;ModG.diag().zeros();
  arma::uvec q=find(ModG==1);
  return q;
}

// [[Rcpp::export]]
arma::mat  Transitive(const arma::mat& G,unsigned int K){
  arma::mat G_red=G;
  arma::mat R_c=Reachability(G_red,K);
  arma::mat C_c=ConnectMat(R_c,K);
  for(unsigned int i=0;i<K;i++){
    for(unsigned int j=0;j<K;j++){
      if(G(i,j)==1){
        arma::mat Gtemp=G_red;Gtemp(i,j)=0;
        arma::mat Rtemp=Reachability(Gtemp,K);
        arma::mat Ctemp=ConnectMat(Rtemp,K);
        if(arma::all(arma::vectorise(Ctemp==C_c))){
          G_red=Gtemp;
        }
      }
    }
  }
  return(G_red);
}

// [[Rcpp::export]]
arma::vec PossibleClass(const arma::mat& R,unsigned int K){//generate all possible arttribute profiles
  arma::vec ClassInd;
  arma::vec vv=bijectionvector(K);
  unsigned int m=K;
  arma::mat AllC=R;//AllC.insert_cols(0,1);
  ClassInd=AllC.t()*vv;
  for(unsigned int j=0;j<K;j++){
    for(unsigned int i=j+1;i<m;i++){
      arma::vec temp=Boolean(AllC.col(j)+AllC.col(i));
      double c_temp=sum(temp%vv);
      if(any((ClassInd==c_temp))){
         continue;
      }else{
        AllC.insert_cols(AllC.n_cols,temp);
        ClassInd=AllC.t()*vv;
        m=AllC.n_cols;
      }
    }
  }
  ClassInd.insert_rows(m,1);
  return ClassInd;
}


// [[Rcpp::export]]
Rcpp::List Collapsed(const arma::mat& R,unsigned int K){//generate all possible arttribute profiles
  arma::vec vv=bijectionvector(K);
  unsigned int NC=pow(2,K);
  arma::vec NewDelta(NC);NewDelta.zeros();
  arma::vec NewInd(NC);NewInd.zeros();
  for(unsigned int c=0;c<NC;c++){
   double cc=c;
   arma::vec aa=inv_bijectionvector(K,cc);
   arma::vec newaa=aa;
   for(unsigned int j=0;j<K;j++){
     if(aa(j)==1){
       arma::vec prereq=R.col(j);
       newaa=Booleanvec(prereq+newaa);
     }
   }
   double ncc=sum(newaa%vv);
   NewInd(c)=ncc;
   NewDelta(ncc)++;
  }
  return Rcpp::List::create(Rcpp::Named("NewInd",NewInd),
                            Rcpp::Named("NewDelta",NewDelta)
  );
}


// [[Rcpp::export]]
double Gammaratio(double L1, double L2){ //compute Gamma ratio
  double ratio=1.0;
  if(L1>L2){
     for(int l=(L1-1);l>L2;l--){
       ratio=ratio*l;
     } 
    }
  if(L1<L2){
      for(int l=(L2-1);l>L1;l--){
        ratio=ratio/l;
      } 
    }
  if(ratio==arma::datum::inf) ratio=MAX;
  return ratio;
}


// [[Rcpp::export]]
double p_ratio_multi(unsigned int K,
                     const arma::vec& Countold, const arma::vec& Countnext,
                     const arma::vec& Deltaold, const arma::vec& Deltanext){
  double ratio=1.0;
  arma::uvec indnext=sort_index(Countnext,"descend");
  arma::vec cnext=Countnext(indnext);arma::vec dnext=Countnext(indnext);
  arma::uvec indold=sort_index(Countold,"descend");
  arma::vec cold=Countold(indold);arma::vec dold=Countold(indold);
 
 /* for(unsigned int c=0;c<NC;c++){
     double ratio1=Gammaratio(cnext(c)+dnext(c),cold(c)+dold(c));
     double temp_ratio=ratio*ratio1;
     if(temp_ratio==temp_ratio& temp_ratio>0){
       ratio=temp_ratio;
     }else{
       break;
       }
  }
   if(ratio==arma::datum::inf) ratio=MAX;
   return ratio;
  */

  double log1=sum((cold+dold+1)%log(cold+dold+1));
  double log2=sum((cnext+dnext+1)%log(cnext+dnext+1));
 // std::cout << "log1:"<<log1<< std::endl;
 // std::cout << "log2:"<<log2<< std::endl;
  ratio=exp(log2-log1); 
 // std::cout << "gammaratio:"<<ratio<< std::endl;
  return ratio;
}


// [[Rcpp::export]]
double p_pi_ratio(unsigned int K, const arma::vec& pis,
                     const arma::mat& Rold, const arma::mat& Rnext){
  double ratio=1.0;
  unsigned int NC=pow(2,K);
  Rcpp::List OLD=Collapsed(Rold,K);
  Rcpp::List NEXT=Collapsed(Rnext,K);
  arma::vec Deltaold=Rcpp::as<arma::vec>(OLD[1]);
  arma::vec Deltanext=Rcpp::as<arma::vec>(NEXT[1]);
  arma::vec INDOLD=Rcpp::as<arma::vec>(OLD[0]);
  arma::vec INDnext=Rcpp::as<arma::vec>(NEXT[0]);
  
  double log1=0;double log2=0;
  arma::vec pisnext(NC);pisnext.zeros();
  for(unsigned int c=0;c<NC;c++){
        double cc=INDnext(c);
        pisnext(cc)=pisnext(cc)+pis(cc);
  }
  log1=sum((Deltaold+1)%log(pis));
  log2=sum((Deltanext+1)%log(pisnext));
  // std::cout << "log1:"<<log1<< std::endl;
  // std::cout << "log2:"<<log2<< std::endl;
  ratio=exp(log2-log1); 
  // std::cout << "gammaratio:"<<ratio<< std::endl;
  return ratio;
  }


// [[Rcpp::export]]
double p_alpha_ratio(const arma::mat& ETA, const arma::mat& Y,
                       const arma::vec& ss,const arma::vec& gs,
                       const arma::mat& Rold,
                       const arma::mat& Rnext,
                       const arma::vec& Class,
                       const arma::vec& pis,
                       unsigned int K,unsigned int N){
  double ratio=1.0;double t2=1.0;
  unsigned int NC=pow(2,K);
  arma::vec countold(NC);
  arma::vec countnext(NC);
  countold.zeros();countnext.zeros();
  Rcpp::List OLD=Collapsed(Rold,K);
  Rcpp::List NEXT=Collapsed(Rnext,K);
  arma::vec Deltaold=Rcpp::as<arma::vec>(OLD[1]);
  arma::vec Deltanext=Rcpp::as<arma::vec>(NEXT[1]);
  arma::vec INDOLD=Rcpp::as<arma::vec>(OLD[0]);
  arma::vec INDnext=Rcpp::as<arma::vec>(NEXT[0]);
  for(unsigned int i=0;i<N;i++){
    arma::vec Yi = (Y.row(i)).t();
    double ci=Class(i);
    double ciold=ci;
    double cinext=INDnext(ci);
    arma::vec PSold(NC);PSold.zeros();
    arma::vec PSnext(NC);PSnext.zeros();
    for(unsigned int cc=0;cc<NC;cc++){
      double ci_old=INDOLD(cc); 
      double ci_next=INDnext(cc);
      arma::vec ETA_i1=ETA.col(ci_old);
      arma::vec ETA_i2 = ETA.col(ci_next);
      PSold(cc)=pYi(ETA_i1,Yi,ss,gs);
      PSnext(cc)=pYi(ETA_i2,Yi,ss,gs);
    }
    t2=sum(PSnext%pis)/sum(PSold%pis);
    ratio=ratio*t2;
 countold(ciold)++;
 countnext(cinext)++;
}
 // std::cout << "Pratio:"<<ratio<< std::endl;
  double ratio1=1.0;
  ratio1=p_pi_ratio(K,pis,Rold,Rnext);
// return (ratio1*ratio);
return ratio;
}




// [[Rcpp::export]]
arma::vec  update_G(unsigned int N,const int K,unsigned int J,
                       const arma::mat& ETA,const arma::mat& Y, 
                       const arma::vec& CLASS,const arma::vec&ss,const arma::vec& gs,
                       const arma::vec& PIsmax,
                       arma::mat& G, arma::mat& R, 
                       double p1,double p2){ 
  //add (1) with p1 or delete (2) with p2
  double move=0.0;
  double add=0.0;
  double ratio=0.0;double ratio2=0.0;
  double r1=0.0;
  double t1=0.0;
  double t2=0.0;
  double ku= R::runif(0.0,1.0);
  if(ku<p1){
      add=1.0;//add
    }else if(ku>(1-p2)){
      add=2.0;//remove
    }else{
      add=0.0;
    }
    arma::mat Rold=R;
    arma::mat Cold=ConnectMat(Rold,K);
if(add==1.0){  //Adding a link
      arma::uvec q=check_no_connect(Cold,K);
      if(q.n_elem<1){
        move=0.0;
      }else{
        arma::uvec temp= arma::randi<arma::uvec>(1,arma::distr_param(0,q.n_elem-1));
        unsigned int i=temp(0);
         arma::mat Gnext=G;
         Gnext(q(i))=1;
         arma::mat Rnext=Reachability(Gnext,K);
         arma::vec Classnext=PossibleClass(Rnext,K);
         arma::mat Cnext=ConnectMat(Rnext,K);
         //make Gnext a transitive reduction:
         Gnext=Transitive(Gnext,K);
         //MH ratio
         r1=p_alpha_ratio(ETA,Y,ss,gs,Rold,Rnext,CLASS,PIsmax,K,N);
    //    std::cout << "add:"<<r1<< std::endl;
         arma::uvec qtemp=check_link(Gnext,K);
         double n2=qtemp.n_elem;
         double n1=q.n_elem;
         t1=1/n1;
         t2=1/n2;
         ratio=r1*p2*t2/(p1*t1);
         if(ratio==arma::datum::inf) ratio=MAX;
         ratio2=fmin(1,ratio);
         double u=R::runif(0,1);
         if(u<ratio2){
           G=Gnext;
           R=Rnext;
           move=1.0;
         }else{
           move=0.0;
           }
       }
   }else if(add==2.0){ //Deleting a link
   arma::uvec q=check_link(G,K);
     if(q.n_elem<1){
       move=0.0;
     }else{
       arma::uvec temp= arma::randi<arma::uvec>(1,arma::distr_param(0,q.n_elem-1));
       unsigned int l=temp(0);
       arma::mat Gnext=G;
       Gnext(q(l))=0;
       arma::mat Rnext=Reachability(Gnext,K);
       arma::vec Classnext=PossibleClass(Rnext,K);
       //MH ratio
       r1=p_alpha_ratio(ETA,Y,ss,gs,Rold,Rnext,CLASS,PIsmax,K,N);
     //  std::cout << "delete:"<<r1<< std::endl;
       arma::mat Cnext=ConnectMat(Rnext,K);
       arma::uvec qtemp=check_no_connect(Cnext,K);
       double n1=qtemp.n_elem;
       double n2=q.n_elem;
       t1=1/n1;
       t2=1/n2;
       ratio=r1*p1*t1/(p2*t2);ratio2=fmin(1,ratio);
       double u=R::runif(0,1);
       if(u<ratio2){
         G=Gnext;
         R=Rnext;
         move=2.0;
       }else{
         move=0.0;
         }
     }
  }
 arma::vec track(3);track(0)=add;track(1)=move;
 track(2)=ratio;
 return(track);
}


// [[Rcpp::export]]
void  parm_update_K(unsigned int N,unsigned int J,unsigned int K,
                    const arma::mat& ETA, const arma::mat& Y, const arma::mat& R,
                    arma::vec& gs,arma::vec& ss,arma::vec& CLASS, 
                    arma::vec& PIsmax){
  unsigned int nClass=pow(2,K);
  arma::vec pY(nClass);pY.zeros();
  arma::cube ab_tilde = arma::zeros<arma::cube>(J,2,2);
  Rcpp::List Match=Collapsed(R,K);
  arma::vec Delta=Rcpp::as<arma::vec>(Match[1]);
  arma::vec IND=Rcpp::as<arma::vec>(Match[0]);
  arma::vec Count(nClass);Count.zeros();
  //update alpha
  double class_i=0;
    for(unsigned int i=0;i<N;i++){
        arma::vec Yi = (Y.row(i)).t();
        for(unsigned int c=0;c<nClass;c++){
            double cc=IND(c);
            arma::vec ETA_i = ETA.col(cc);
            pY(cc) = pYi(ETA_i,Yi,ss,gs);
        }
        arma::vec numerator = pY%(Count+Delta);
        arma::vec PS = numerator/arma::sum(numerator);
        class_i = rmultinomial(PS);
        CLASS(i) = class_i;
        Count(class_i)++;
        //update guess and slip full conditional beta parms
        arma::vec ETA_c = ETA.col(class_i);
        for(unsigned int j=0;j<J;j++){
            ab_tilde(j,ETA_c(j),Yi(j)) += 1.;
        }
    }
    arma::vec deltatilde = Count;
    arma::vec pis = rDirichlet(deltatilde+Delta);
    PIsmax=pis;
  //update guess and slip probabilities
  for(unsigned int j=0;j<J;j++){
    double us=R::runif(0,1);
    double ug=R::runif(0,1);
    double sold = ss(j);
    //draw g conditoned upon s_t-1
    double ab_g1 = ab_tilde(j,0,1);
    double ab_g0 = ab_tilde(j,0,0);
    double pg = R::pbeta(1.0-sold,ab_g1+1.,ab_g0+1.,1,0);
    double gnew = R::qbeta(ug*pg,ab_g1+1.,ab_g0+1.,1,0);
    //draw s conditoned upon g
    double ab_s1 = ab_tilde(j,1,1);
    double ab_s0 = ab_tilde(j,1,0);
    double ps = R::pbeta(1.0-gnew,ab_s0+1.,ab_s1+1.,1,0);
    double snew = R::qbeta(us*ps,ab_s0+1.,ab_s1+1.,1,0);
    gs(j) = gnew;
    ss(j) = snew;
  }
}



// [[Rcpp::export]]
Rcpp::List dina_AHM(const arma::mat& Y, const arma::mat& Q,const unsigned int K, 
                    double p1,double p2,
                    unsigned int burnin=10000, unsigned int chain_length=15000){
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  unsigned int chain_m_burn = chain_length-burnin;
  unsigned int tmburn;
  unsigned int nClass = pow(2,K);
  //Savinging output
  arma::mat PIs(nClass,chain_m_burn);PIs.zeros();
  arma::mat SS(J,chain_m_burn);
  arma::mat GS(J,chain_m_burn);
  arma::mat Class(N,chain_m_burn);
  arma::cube Graphs(K,K,chain_m_burn);
  arma::cube RS(K,K,chain_m_burn);
  arma::mat adds(3,chain_length);
  arma::vec m2LL(chain_m_burn);
  
    //initialize theta,pis 
  arma::vec ss = arma::randu<arma::vec>(J);
  arma::vec gs = (arma::ones<arma::vec>(J) - ss)%arma::randu<arma::vec>(J);

  //initialize G
  arma::mat G(K,K);G.zeros();
  arma::mat R=Reachability(G,K);
  arma::mat C=ConnectMat(R,K);
  arma::vec ClassInd=PossibleClass(R,K);
  arma::vec CLASS=arma::randi<arma::vec>(N,arma::distr_param(0,nClass-1));
  arma::vec delta0 = arma::ones<arma::vec>(nClass);

  arma::vec Countold(nClass);Countold.zeros();
  for(unsigned i=0;i<N;i++){
    double c=CLASS(i);
    Countold(c)++;
  }
  arma::vec PIsmax = rDirichlet(delta0);  
      
  arma::mat ETA = ETAmat(K,J,Q);
  arma::vec vv = bijectionvector(K);
  
  //Start Markov chain
  //unsigned int unq=0; //start with 0 mode
  for(unsigned int t = 0; t < chain_length; t++){
    //updata G,R
   adds.col(t)=update_G(N,K,J,ETA,Y,CLASS,ss,gs,PIsmax,G,R,p1,p2);
    //update parm via pointer
   parm_update_K(N,J,K,ETA,Y,R,gs,ss,CLASS,PIsmax);
    
    //unsigned int flag=0; //not mode
    if(t>burnin-1){
      double m2ll=dina_m2LL(N,J,K,Y,PIsmax,Q,ss,gs);
      tmburn = t-burnin;
      //update parameter value via pointer. save classes and PIs
      SS.col(tmburn)       = ss;
      GS.col(tmburn)       = gs;
      Class.col(tmburn)    =CLASS;
      Graphs.slice(tmburn) =G;
      RS.slice(tmburn)     =R;
      PIs.col(tmburn)      = PIsmax;
      m2LL(tmburn)         = m2ll;
    }
  if (t%1000 == 0) {
     // std::cout << t << std::endl;
  }
  }
  return Rcpp::List::create(Rcpp::Named("GS",GS),
                            Rcpp::Named("SS",SS),
                            Rcpp::Named("PIs",PIs),
                            Rcpp::Named("CLASS",Class),
                            Rcpp::Named("RS",RS),
                            Rcpp::Named("Graphs",Graphs),
                            Rcpp::Named("Trace",adds),
                            Rcpp::Named("m2LL",m2LL)
  );
}



// [[Rcpp::export]]
Rcpp::List dina_AHM_check(const arma::mat& Y, const arma::mat& Q,const unsigned int K, 
                    double p1,double p2,
                    const arma::vec CLASS,
                    const arma::vec ss,
                    const arma::vec gs,const arma::vec PIsmax,
                    unsigned int burnin=10000, unsigned int chain_length=15000){
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  unsigned int chain_m_burn = chain_length-burnin;
  unsigned int tmburn;
  unsigned int nClass = pow(2,K);
  //Savinging output
 // arma::mat PIs(nClass,chain_m_burn);PIs.zeros();
//  arma::mat SS(J,chain_m_burn);
//  arma::mat GS(J,chain_m_burn);
//  arma::mat Class(N,chain_m_burn);
  arma::cube Graphs(K,K,chain_m_burn);
  arma::cube RS(K,K,chain_m_burn);
  arma::mat adds(3,chain_length);
  arma::vec m2LL(chain_m_burn);
  
  //initialize theta,pis 
//  arma::vec ss = arma::randu<arma::vec>(J);
//  arma::vec gs = (arma::ones<arma::vec>(J) - ss)%arma::randu<arma::vec>(J);
  
  //initialize G
  arma::mat G(K,K);G.zeros();
  arma::mat R=Reachability(G,K);
  arma::mat C=ConnectMat(R,K);
//  arma::vec ClassInd=PossibleClass(R,K);
//  arma::vec CLASS=arma::randi<arma::vec>(N,arma::distr_param(0,nClass-1));
  arma::vec delta0 = arma::ones<arma::vec>(nClass);
  
  arma::vec Countold(nClass);Countold.zeros();
  for(unsigned i=0;i<N;i++){
    double c=CLASS(i);
    Countold(c)++;
  }
 // arma::vec PIsmax = rDirichlet(delta0);  
  
  arma::mat ETA = ETAmat(K,J,Q);
  arma::vec vv = bijectionvector(K);
  
  //Start Markov chain
  //unsigned int unq=0; //start with 0 mode
  for(unsigned int t = 0; t < chain_length; t++){
    //updata G,R
    adds.col(t)=update_G(N,K,J,ETA,Y,CLASS,ss,gs,PIsmax,G,R,p1,p2);
    //update parm via pointer
 //   parm_update_K(N,J,K,ETA,Y,R,gs,ss,CLASS,PIsmax);
    
    //unsigned int flag=0; //not mode
    if(t>burnin-1){
      double m2ll=dina_m2LL(N,J,K,Y,PIsmax,Q,ss,gs);
      tmburn = t-burnin;
      //update parameter value via pointer. save classes and PIs
   //   SS.col(tmburn)       = ss;
    //  GS.col(tmburn)       = gs;
    //  Class.col(tmburn)    =CLASS;
      Graphs.slice(tmburn) =G;
      RS.slice(tmburn)     =R;
      //PIs.col(tmburn)      = PIsmax;
      m2LL(tmburn)         = m2ll;
    }
    if (t%1000 == 0) {
   //   std::cout << t << std::endl;
    }
  }
  return Rcpp::List::create(//Rcpp::Named("GS",GS),
                           // Rcpp::Named("SS",SS),
                          //  Rcpp::Named("PIs",PIs),
                          //  Rcpp::Named("CLASS",Class),
                            Rcpp::Named("RS",RS),
                            Rcpp::Named("Graphs",Graphs),
                            Rcpp::Named("Trace",adds),
                            Rcpp::Named("m2LL",m2LL)
  );
}


