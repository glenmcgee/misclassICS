#include <Rcpp.h>
using namespace Rcpp;

// May 27, 2018 EDIT:
// sigma and gamma are vectors of length K (because for random slopes it differs by cluster)


// expit function for use in C code
double expit(double x){
  double out = exp(x)/(1+exp(x));
  return out;
}

// expit function for use in C code
long double max(double x, double y){
  long double out;
  if(x>y){
    out = x;
  } else {
    out = y;
  }
  return out;
}

// expit function for use in C code
NumericMatrix square_vect(NumericVector x){
  int p = x.size(); // dimensions
  NumericMatrix out(p,p);
  int i, j;

  for(i = 0; i < p; i++){ // loop over rows

    for(j = 0; j < p; j++){ // loop over cols

      out(i,j) = x[i]*x[j];

    }

  }

  return out;
}



// [[Rcpp::export]]
NumericVector compute_Omega(NumericVector gamma,
                            NumericVector sigma,
                            NumericVector zeta,
                            NumericVector GH_w, // these include scaling factors
                            NumericVector GH_z) {

  int maxK = zeta.size(); // number of clusters
  int nquad = GH_w.size(); // number of quadrature points
  NumericVector Omega(maxK);
  int iter, converged, q;
  double zetak, Omegak, g_Omegak, g_zetak, step;

  for(int k = 0; k < maxK; k++) { // loop over clusters

    zetak = zeta[k];
    Omegak = zetak; // starting value

    g_zetak = exp(zetak);  // marginal value

    iter = 0;
    converged = 0;

    do{

      g_Omegak = 0;

      for(q = 0; q < nquad; q++){ // loop over GH zeros

        // this is the integral side of the function to be solved
        // also equivalent to A1k
        g_Omegak += GH_w[q]*exp(Omegak+gamma[k]*sigma[k]*GH_z[q]);

      }

      //NR step is function/gradient
      step = (g_Omegak-g_zetak)/(g_Omegak+1e-323);   // 1e-323 is for singularity
      if (step < -0.5) step = -0.5;               // in case of huge steps
      if (step > 0.5) step = 0.5;

      Omegak = Omegak - step; // Newton Raphson

      if( fabs(step) < 1e-7 ) converged=1;
      iter++;

    } while(iter<200 && !converged );

    Omega[k] = Omegak;
  }
  return Omega;
}


// [[Rcpp::export]]
NumericVector compute_Delta(NumericVector sigma,
                            NumericVector eta,
                            NumericVector GH_w, // these include scaling factors
                            NumericVector GH_z) {

  int N = eta.size(); // number of outcomes overall
  int nquad = GH_w.size(); // number of quadrature points
  NumericVector Delta(N);
  int iter, converged, q;
  double etaki, Deltaki, h_Deltaki, h_etaki, A2ki, step;


  for(int ki = 0; ki < N; ki++) { // loop over all outcomes

    etaki = eta[ki];
    Deltaki = etaki;

    h_etaki = expit(etaki); // marginal side

    iter = 0;
    converged = 0;

    do{

      h_Deltaki = 0;
      A2ki = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros

        // this is the integral side of the function to be solved
        h_Deltaki += GH_w[q]*expit(Deltaki+sigma[ki]*GH_z[q]);
        // gradient
        A2ki += GH_w[q]*(expit(Deltaki+sigma[ki]*GH_z[q]))*(1-expit(Deltaki+sigma[ki]*GH_z[q]));

      }

      //NR step is function/gradient
      step = (h_Deltaki-h_etaki)/(A2ki+1e-323);   // 1e-323 is for singularity
      if (step< -0.5) step = -0.5;               // in case of huge steps
      if (step> 0.5) step = 0.5;

      Deltaki = Deltaki - step; // Newton Raphson

      if( fabs(step) < 1e-7 ) converged=1;
      iter++;

    } while(iter<200 && !converged );

    Delta[ki] = Deltaki;
  }
  return Delta;
}


// [[Rcpp::export]]
double compute_logLik_MISclassICS(NumericVector gamma,
                                  NumericVector gamma0,
                                  NumericVector gamma1,
                                  NumericVector sigma,
                                  NumericVector sigma0,
                                  NumericVector sigma1,
                                  NumericVector Omega,
                                  NumericVector Omega0,
                                  NumericVector Omega1,
                                  NumericVector Delta,
                                  NumericVector Delta0,
                                  NumericVector Delta1,
                                  NumericVector etaW,
                                  NumericVector etaW0,
                                  NumericVector etaW1,
                                  NumericVector etaX,
                                  NumericVector Nk,
                                  NumericVector Nk_actual,
                                  NumericVector Yki,
                                  NumericVector Wk,
                                  NumericVector Xk,
                                  NumericVector IDk,
                                  LogicalVector VALIDk,
                                  int minNk,
                                  NumericVector GH_w,
                                  NumericVector GH_z) {

  int maxK = max(IDk); // number of clusters
  int nquad = GH_w.size(); // number of quadrature points
  int cumNk = 0; // to only subset outcomes from cluster k
  long double logL = 0; // log Likelihood (objective function to be maximized)
  int k, q, i;
  long double L_k, L_YkNk, L_Nk, L_Yk, L_YkNk0, L_Nk0, L_Yk0, L_YkNk1, L_Nk1, L_Yk1;
  long double L_Wk, L_Wk0, L_Wk1;
  long double L_Xk, P_Xk1;
  NumericVector Nk_out = Nk - minNk; // adjust for Poisson + minimum
  NumericVector logNk_factorial = lgamma(Nk_out+1);

  for(k = 0; k < maxK; k++) { // loop over clusters
    L_YkNk = 0; // kth contribution to the (non-log) likelihood
    L_YkNk0 = 0;
    L_YkNk1 = 0;

    for(q = 0; q < nquad; q++){ // loop over GH zeros
      // its really logLNk
      L_Nk =  Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q]) -exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k];   // contribution of Nk to likelihood
      L_Nk0 =  Nk_out[k]*(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]) -exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q])-logNk_factorial[k];
      L_Nk1 =  Nk_out[k]*(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]) -exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q])-logNk_factorial[k];
      // contribution of Yk to likelihood
      L_Yk = 1;
      L_Yk0 = 1;
      L_Yk1 = 1;
      for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
        L_Yk *=  exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q]));
        L_Yk0 *=  exp(Yki[i]*(Delta0[i]+sigma0[k]*GH_z[q]))/(1+exp(Delta0[i]+sigma0[k]*GH_z[q]));
        L_Yk1 *=  exp(Yki[i]*(Delta1[i]+sigma1[k]*GH_z[q]))/(1+exp(Delta1[i]+sigma1[k]*GH_z[q]));
      }
      L_YkNk +=  exp(log(GH_w[q])+L_Nk+log(L_Yk));
      L_YkNk0 +=  exp(log(GH_w[q])+L_Nk0+log(L_Yk0));
      L_YkNk1 +=  exp(log(GH_w[q])+L_Nk1+log(L_Yk1));
    }
    L_Wk =  exp(Wk[k]*(etaW[k]))/(1+exp(etaW[k])); // misclassification likelihood
    L_Wk0 =  exp(Wk[k]*(etaW0[k]))/(1+exp(etaW0[k]));
    L_Wk1 =  exp(Wk[k]*(etaW1[k]))/(1+exp(etaW1[k]));
    // single exposure likelihood
    L_Xk =  exp(Xk[k]*(etaX[k]))/(1+exp(etaX[k]));
    P_Xk1 =  exp(1*(etaX[k]))/(1+exp(etaX[k]));

    if(VALIDk[k] == TRUE){
      L_k = L_Wk*L_YkNk*L_Xk;
    } else {
      L_k = L_Wk0*L_YkNk0*(1-P_Xk1)+L_Wk1*L_YkNk1*P_Xk1;
    }

    cumNk += Nk_actual[k];

    logL += log(L_k+1e-323);
  }

  return logL;
}


// note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// [[Rcpp::export]]
NumericVector compute_score_MISclassICS(NumericVector gamma,
                                        NumericVector gamma0,
                                        NumericVector gamma1,
                                        NumericVector sigma,
                                        NumericVector sigma0,
                                        NumericVector sigma1,
                                        NumericVector Omega,
                                        NumericVector Omega0,
                                        NumericVector Omega1,
                                        NumericVector Delta,
                                        NumericVector Delta0,
                                        NumericVector Delta1,
                                        NumericVector zeta,
                                        NumericVector zeta0,
                                        NumericVector zeta1,
                                        NumericVector eta,
                                        NumericVector eta0,
                                        NumericVector eta1,
                                        NumericVector etaW,
                                        NumericVector etaW0,
                                        NumericVector etaW1,
                                        NumericVector etaX,
                                        NumericVector Nk,
                                        NumericVector Nk_actual,
                                        NumericMatrix Zk,
                                        NumericMatrix Zk0,
                                        NumericMatrix Zk1,
                                        NumericVector Yki,
                                        NumericMatrix Xki,
                                        NumericMatrix Xki0,
                                        NumericMatrix Xki1,
                                        NumericVector Wk,
                                        NumericMatrix Zwk,
                                        NumericMatrix Zwk0,
                                        NumericMatrix Zwk1,
                                        NumericVector Xk,
                                        NumericMatrix Zxk,
                                        NumericVector IDk,
                                        LogicalVector VALIDk,
                                        int minNk,
                                        bool condSize,
                                        bool condOut,
                                        NumericVector GH_w,
                                        NumericVector GH_z,
                                        int len_alpha,
                                        int len_beta,
                                        int len_epsW,
                                        int len_epsX) {

  int maxK = Omega.size(); // number of clusters
  int N = Delta.size(); // number of obs
  int nquad = GH_w.size(); // number of quadrature points
  int cumNk = 0; // to only subset outcomes from cluster k
  int k, q, i, j;
  NumericVector Nk_out = Nk - minNk;
  NumericVector logNk_factorial = lgamma(Nk_out+1);
  NumericVector score(len_alpha+len_beta+1+1+len_epsW+len_epsX);
  long double L_k;
  long double L_YkNk, L_Nk, L_Yk, A1k, B1k, L_YkNk0, L_Nk0, L_Yk0, A1k0, B1k0, L_YkNk1, L_Nk1, L_Yk1, A1k1, B1k1;
  long double L_Wk, L_Xk, L_Wk0, L_Xk0, L_Wk1, L_Xk1;
  NumericVector A2ki(N), B2ki(N), A2ki0(N), B2ki0(N), A2ki1(N), B2ki1(N);
  long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma, dOmegak_dzetak0, dOmegak_dgamma0, dOmegak_dsigma0, dOmegak_dzetak1, dOmegak_dgamma1, dOmegak_dsigma1;
  NumericVector dOmegak_dalpha(len_alpha), dOmegak_dalpha0(len_alpha), dOmegak_dalpha1(len_alpha);
  NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N), ddeltaki_detaki0(N), ddeltaki_sigma0(N), ddeltaki_detaki1(N), ddeltaki_sigma1(N);
  NumericMatrix ddeltaki_dbeta(N,len_beta), ddeltaki_dbeta0(N,len_beta), ddeltaki_dbeta1(N,len_beta);
  long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma, dLNk_dgamma0, dLNk_dsigma0, dLYk_dsigma0, dLk_dgamma0, dLk_dsigma0, dLNk_dgamma1, dLNk_dsigma1, dLYk_dsigma1, dLk_dgamma1, dLk_dsigma1;
  NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha), dLNk_dalpha0(len_alpha), dLk_dalpha0(len_alpha), dLNk_dalpha1(len_alpha), dLk_dalpha1(len_alpha);
  NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta), dLYk_dbeta0(len_beta), dLk_dbeta0(len_beta), dLYk_dbeta1(len_beta), dLk_dbeta1(len_beta);
  NumericVector dLk_depsW(len_epsW), dLk_depsX(len_epsX), dLk_depsW0(len_epsW), dLk_depsX0(len_epsX), dLk_depsW1(len_epsW), dLk_depsX1(len_epsX);

  for(k = 0; k < maxK; k++) { // loop over clusters


    if(VALIDk[k] == TRUE){

      A1k = 0;
      B1k = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
        B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];

        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
          B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
        }
      }

      // at this point we have A1 B1 A2 B2, eta, zeta

      // derivative chunks for size model
      if (condSize == FALSE) {
        dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
        dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
        dOmegak_dzetak = exp(zeta[k])/(A1k+1e-323);
      } else{  // much simpler if conditional model is used
        dOmegak_dgamma = 0;
        dOmegak_dsigma = 0;
        dOmegak_dzetak = 1;
      }
      for(j = (0); j < (len_alpha); j++){
        dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
      }

      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_k = 0;
      L_YkNk = 0;
      dLk_dalpha = 0*dLk_dalpha;
      dLk_dgamma = 0;
      dLk_dbeta = 0*dLk_dbeta;
      dLk_dsigma = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        // likelihood and derivative contributions form size
        L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
        for(j = (0); j < (len_alpha); j++){
          dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
        }
        dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
        dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));


        // likelihood and derivative contributions from outcome
        L_Yk = 1; // contribution of Yk to likelihood
        dLYk_dsigma = 0;
        dLYk_dbeta = 0*dLYk_dbeta;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          L_Yk *=  exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta[j] += ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
          }
          dLYk_dsigma += (ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
        }
        dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required

        L_YkNk += GH_w[q]*L_Nk*L_Yk;
        // compute main derivs here
        for(j = (0); j < (len_alpha); j++){
          dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
        }
        dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
        }
        dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma);
      }

      L_Wk = exp(Wk[k]*(etaW[k]))/(1+exp(etaW[k])); // misclassification likelihood;
      L_Xk = exp(Xk[k]*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW[j] = L_Wk*Zwk(k,j)*(Wk[k]-expit(etaW[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX[j] = L_Xk*Zxk(k,j)*(Xk[k]-expit(etaX[k]));
      }

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
        score[j] += dLk_dalpha[j]/(L_YkNk+1e-323);
      }
      score[len_alpha] += dLk_dgamma/(L_YkNk+1e-323);
      for(j = 0; j < (len_beta); j++){
        score[j+len_alpha+1] += dLk_dbeta[j]/(L_YkNk+1e-323);
      }
      score[len_alpha+1+len_beta] += dLk_dsigma/(L_YkNk+1e-323);
      for(j = 0; j < (len_epsW); j++){
        score[j+len_alpha+1+len_beta+1] += dLk_depsW[j]/(L_Wk+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_alpha+1+len_beta+1+len_epsW] += dLk_depsX[j]/(L_Xk+1e-323);
      }
    }else{ // main sample

      A1k0 = 0;
      A1k1 = 0;
      B1k0 = 0;
      B1k1 = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        A1k0 += GH_w[q]*exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]);
        A1k1 += GH_w[q]*exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]);
        B1k0 += GH_w[q]*exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q])*GH_z[q];
        B1k1 += GH_w[q]*exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q])*GH_z[q];

        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          A2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          A2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          B2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]))*GH_z[q];
          B2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]))*GH_z[q];
        }
      }

      // at this point we have A1 B1 A2 B2, eta, zeta

      // derivative chunks for size model
      if (condSize == FALSE) {
        dOmegak_dgamma0 = -sigma0[k]*B1k0/(A1k0+1e-323);
        dOmegak_dgamma1 = -sigma1[k]*B1k1/(A1k1+1e-323);
        dOmegak_dsigma0 = -gamma0[k]*B1k0/(A1k0+1e-323);
        dOmegak_dsigma1 = -gamma1[k]*B1k1/(A1k1+1e-323);
        dOmegak_dzetak0 = exp(zeta0[k])/(A1k0+1e-323);
        dOmegak_dzetak1 = exp(zeta1[k])/(A1k1+1e-323);
      } else{  // much simpler if conditional model is used
        dOmegak_dgamma0 = 0;
        dOmegak_dgamma1 = 0;
        dOmegak_dsigma0 = 0;
        dOmegak_dsigma1 = 0;
        dOmegak_dzetak0 = 1;
        dOmegak_dzetak1 = 1;
      }
      for(j = (0); j < (len_alpha); j++){
        dOmegak_dalpha0[j] =  dOmegak_dzetak0*Zk0(k,j);
        dOmegak_dalpha1[j] =  dOmegak_dzetak1*Zk1(k,j);
      }

      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = expit(eta0[i])*(1-expit(eta0[i]))/(A2ki0[i]+1e-323);
          ddeltaki_detaki1[i] = expit(eta1[i])*(1-expit(eta1[i]))/(A2ki1[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = -B2ki0[i]/(A2ki0[i]+1e-323);
          ddeltaki_sigma1[i] = -B2ki1[i]/(A2ki1[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = 1;
          ddeltaki_detaki1[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = 0;
          ddeltaki_sigma1[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_k = 0;
      L_YkNk0 = 0;
      L_YkNk1 = 0;
      dLk_dalpha0 = 0*dLk_dalpha0;
      dLk_dalpha1 = 0*dLk_dalpha1;
      dLk_dgamma0 = 0;
      dLk_dgamma1 = 0;
      dLk_dbeta0 = 0*dLk_dbeta0;
      dLk_dbeta1 = 0*dLk_dbeta1;
      dLk_dsigma0 = 0;
      dLk_dsigma1 = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        // likelihood and derivative contributions form size
        L_Nk0 = exp(Nk_out[k]*(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q])-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
        L_Nk1 = exp(Nk_out[k]*(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q])-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
        for(j = (0); j < (len_alpha); j++){
          dLNk_dalpha0[j] =  L_Nk0*dOmegak_dalpha0[j]*(Nk_out[k]-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]));
          dLNk_dalpha1[j] =  L_Nk1*dOmegak_dalpha1[j]*(Nk_out[k]-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]));
        }
        dLNk_dgamma0 = L_Nk0*(dOmegak_dgamma0+sigma0[k]*GH_z[q])*(Nk_out[k]-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]));
        dLNk_dgamma1 = L_Nk1*(dOmegak_dgamma1+sigma1[k]*GH_z[q])*(Nk_out[k]-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]));
        dLNk_dsigma0 = L_Nk0*(dOmegak_dsigma0+gamma0[k]*GH_z[q])*(Nk_out[k]-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]));
        dLNk_dsigma1 = L_Nk1*(dOmegak_dsigma1+gamma1[k]*GH_z[q])*(Nk_out[k]-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]));

        // likelihood and derivative contributions from outcome
        L_Yk0 = 1; // contribution of Yk to likelihood
        L_Yk1 = 1;
        dLYk_dsigma0 = 0;
        dLYk_dsigma1 = 0;
        dLYk_dbeta0 = 0*dLYk_dbeta0;
        dLYk_dbeta1 = 0*dLYk_dbeta1;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          L_Yk0 *=  exp(Yki[i]*(Delta0[i]+sigma0[k]*GH_z[q]))/(1+exp(Delta0[i]+sigma0[k]*GH_z[q]));
          L_Yk1 *=  exp(Yki[i]*(Delta1[i]+sigma1[k]*GH_z[q]))/(1+exp(Delta1[i]+sigma1[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta0[j] += ddeltaki_dbeta0(i,j)*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
            dLYk_dbeta1[j] += ddeltaki_dbeta1(i,j)*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          }
          dLYk_dsigma0 += (ddeltaki_sigma0[i]+GH_z[q])*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          dLYk_dsigma1 += (ddeltaki_sigma1[i]+GH_z[q])*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta0[j] =  L_Yk0*dLYk_dbeta0[j];  // scaling it by L_Yk as required
          dLYk_dbeta1[j] =  L_Yk1*dLYk_dbeta1[j];
        }
        dLYk_dsigma0 = L_Yk0*dLYk_dsigma0; // scaling it by L_Yk as required
        dLYk_dsigma1 = L_Yk1*dLYk_dsigma1;

        L_YkNk0 += GH_w[q]*L_Nk0*L_Yk0;
        L_YkNk1 += GH_w[q]*L_Nk1*L_Yk1;
        // compute main derivs here
        for(j = (0); j < (len_alpha); j++){
          dLk_dalpha0[j] += GH_w[q]*L_Yk0*dLNk_dalpha0[j];
          dLk_dalpha1[j] += GH_w[q]*L_Yk1*dLNk_dalpha1[j];
        }
        dLk_dgamma0 += GH_w[q]*L_Yk0*dLNk_dgamma0;
        dLk_dgamma1 += GH_w[q]*L_Yk1*dLNk_dgamma1;
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta0[j] += GH_w[q]*dLYk_dbeta0[j]*L_Nk0;
          dLk_dbeta1[j] += GH_w[q]*dLYk_dbeta1[j]*L_Nk1;
        }
        dLk_dsigma0 += GH_w[q]*sigma0[k]*(dLYk_dsigma0*L_Nk0+L_Yk0*dLNk_dsigma0);
        dLk_dsigma1 += GH_w[q]*sigma1[k]*(dLYk_dsigma1*L_Nk1+L_Yk1*dLNk_dsigma1);
      }

      L_Wk0 = exp(Wk[k]*(etaW0[k]))/(1+exp(etaW0[k])); // misclassification likelihood;
      L_Wk1 = exp(Wk[k]*(etaW1[k]))/(1+exp(etaW1[k]));
      L_Xk0 = exp(0*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      L_Xk1 = exp(1*(etaX[k]))/(1+exp(etaX[k]));
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW0[j] = L_Wk0*Zwk0(k,j)*(Wk[k]-expit(etaW0[k]));
        dLk_depsW1[j] = L_Wk1*Zwk1(k,j)*(Wk[k]-expit(etaW1[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX0[j] = L_Xk0*Zxk(k,j)*(0-expit(etaX[k]));
        dLk_depsX1[j] = L_Xk1*Zxk(k,j)*(1-expit(etaX[k]));
      }

      // compute denominator
      L_k = L_Wk0*L_YkNk0*L_Xk0+L_Wk1*L_YkNk1*L_Xk1;

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
        score[j] += (L_Wk0*dLk_dalpha0[j]*L_Xk0+L_Wk1*dLk_dalpha1[j]*L_Xk1)/(L_k+1e-323);
      }
      score[len_alpha] += (L_Wk0*dLk_dgamma0*L_Xk0+L_Wk1*dLk_dgamma1*L_Xk1)/(L_k+1e-323);
      for(j = 0; j < (len_beta); j++){
        score[j+len_alpha+1] += (L_Wk0*dLk_dbeta0[j]*L_Xk0+L_Wk1*dLk_dbeta1[j]*L_Xk1)/(L_k+1e-323);
      }
      score[len_alpha+1+len_beta] += (L_Wk0*dLk_dsigma0*L_Xk0+L_Wk1*dLk_dsigma1*L_Xk1)/(L_k+1e-323);
      for(j = 0; j < (len_epsW); j++){
        score[j+len_alpha+1+len_beta+1] += (dLk_depsW0[j]*L_YkNk0*L_Xk0+dLk_depsW1[j]*L_YkNk1*L_Xk1)/(L_k+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_alpha+1+len_beta+1+len_epsW] += (L_Wk0*L_YkNk0*dLk_depsX0[j]+L_Wk1*L_YkNk1*dLk_depsX1[j])/(L_k+1e-323);
      }
    }


    cumNk += Nk_actual[k];
  }

  return score;
}

// note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// [[Rcpp::export]]
NumericVector compute_score_MISclassICS_slopes(NumericVector gamma,
                                               NumericVector gamma0,
                                               NumericVector gamma1,
                                               NumericVector sigma,
                                               NumericVector sigma0,
                                               NumericVector sigma1,
                                               NumericVector Omega,
                                               NumericVector Omega0,
                                               NumericVector Omega1,
                                               NumericVector Delta,
                                               NumericVector Delta0,
                                               NumericVector Delta1,
                                               NumericVector zeta,
                                               NumericVector zeta0,
                                               NumericVector zeta1,
                                               NumericVector eta,
                                               NumericVector eta0,
                                               NumericVector eta1,
                                               NumericVector etaW,
                                               NumericVector etaW0,
                                               NumericVector etaW1,
                                               NumericVector etaX,
                                               NumericVector rand_indic,
                                               NumericVector Nk,
                                               NumericVector Nk_actual,
                                               NumericMatrix Zk,
                                               NumericMatrix Zk0,
                                               NumericMatrix Zk1,
                                               NumericVector Yki,
                                               NumericMatrix Xki,
                                               NumericMatrix Xki0,
                                               NumericMatrix Xki1,
                                               NumericVector Wk,
                                               NumericMatrix Zwk,
                                               NumericMatrix Zwk0,
                                               NumericMatrix Zwk1,
                                               NumericVector Xk,
                                               NumericMatrix Zxk,
                                               NumericVector IDk,
                                               LogicalVector VALIDk,
                                               int minNk,
                                               bool condSize,
                                               bool condOut,
                                               NumericVector GH_w,
                                               NumericVector GH_z,
                                               int len_alpha,
                                               int len_beta,
                                               int len_epsW,
                                               int len_epsX) {

  int maxK = Omega.size(); // number of clusters
  int N = Delta.size(); // number of obs
  int nquad = GH_w.size(); // number of quadrature points
  int cumNk = 0; // to only subset outcomes from cluster k
  int k, q, i, j;
  NumericVector Nk_out = Nk - minNk;
  NumericVector logNk_factorial = lgamma(Nk_out+1);
  NumericVector score(len_alpha+len_beta+2+2+len_epsW+len_epsX); //2 gammas and 2 sigmas
  long double  L_k;
  long double L_YkNk, L_Nk, L_Yk, A1k, B1k, L_YkNk0, L_Nk0, L_Yk0, A1k0, B1k0, L_YkNk1, L_Nk1, L_Yk1, A1k1, B1k1;
  long double L_Wk, L_Xk, L_Wk0, L_Xk0, L_Wk1, L_Xk1;
  NumericVector A2ki(N), B2ki(N), A2ki0(N), B2ki0(N), A2ki1(N), B2ki1(N);
  long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma, dOmegak_dzetak0, dOmegak_dgamma0, dOmegak_dsigma0, dOmegak_dzetak1, dOmegak_dgamma1, dOmegak_dsigma1;
  NumericVector dOmegak_dalpha(len_alpha), dOmegak_dalpha0(len_alpha), dOmegak_dalpha1(len_alpha);
  NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N), ddeltaki_detaki0(N), ddeltaki_sigma0(N), ddeltaki_detaki1(N), ddeltaki_sigma1(N);
  NumericMatrix ddeltaki_dbeta(N,len_beta), ddeltaki_dbeta0(N,len_beta), ddeltaki_dbeta1(N,len_beta);
  long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma, dLNk_dgamma0, dLNk_dsigma0, dLYk_dsigma0, dLk_dgamma0, dLk_dsigma0, dLNk_dgamma1, dLNk_dsigma1, dLYk_dsigma1, dLk_dgamma1, dLk_dsigma1;
  NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha), dLNk_dalpha0(len_alpha), dLk_dalpha0(len_alpha), dLNk_dalpha1(len_alpha), dLk_dalpha1(len_alpha);
  NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta), dLYk_dbeta0(len_beta), dLk_dbeta0(len_beta), dLYk_dbeta1(len_beta), dLk_dbeta1(len_beta);
  NumericVector dLk_depsW(len_epsW), dLk_depsX(len_epsX), dLk_depsW0(len_epsW), dLk_depsX0(len_epsX), dLk_depsW1(len_epsW), dLk_depsX1(len_epsX);

  for(k = 0; k < maxK; k++) { // loop over clusters

    if(VALIDk[k] == TRUE){

      A1k = 0;
      B1k = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros

        A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);

        B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];

        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes

          A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));

          B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];

        }

      }

      // at this point we have A1 B1 A2 B2, eta, zeta

      // derivative chunks for size model
      if (condSize == FALSE) {
        dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
        dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
        dOmegak_dzetak = exp(zeta[k])/(A1k+1e-323);
      } else{  // much simpler if conditional model is used
        dOmegak_dgamma = 0;
        dOmegak_dsigma = 0;
        dOmegak_dzetak = 1;
      }
      for(j = (0); j < (len_alpha); j++){
        dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
      }

      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_YkNk = 0;
      dLk_dalpha = 0*dLk_dalpha;
      dLk_dgamma = 0;
      dLk_dbeta = 0*dLk_dbeta;
      dLk_dsigma = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros

        // likelihood and derivative contributions form size
        L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
        for(j = (0); j < (len_alpha); j++){
          dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
        }
        dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
        dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));


        // likelihood and derivative contributions from outcome
        L_Yk = 1; // contribution of Yk to likelihood
        dLYk_dsigma = 0;
        dLYk_dbeta = 0*dLYk_dbeta;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes

          L_Yk *=  exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta[j] += ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
          }
          dLYk_dsigma += (ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
        }
        dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required

        L_YkNk += GH_w[q]*L_Nk*L_Yk;
        // compute main derivs here
        for(j = (0); j < (len_alpha); j++){
          dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
        }
        dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
        }
        dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma);

      }
      L_Wk = exp(Wk[k]*(etaW[k]))/(1+exp(etaW[k])); // misclassification likelihood;
      L_Xk = exp(Xk[k]*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW[j] = L_Wk*Zwk(k,j)*(Wk[k]-expit(etaW[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX[j] = L_Xk*Zxk(k,j)*(Xk[k]-expit(etaX[k]));
      }

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
        score[j] += dLk_dalpha[j]/(L_YkNk+1e-323);
      }
      score[len_alpha] += (1-rand_indic[k])*dLk_dgamma/(L_YkNk+1e-323);
      score[len_alpha+1] += rand_indic[k]*dLk_dgamma/(L_YkNk+1e-323);
      for(j = 0; j < (len_beta); j++){
        score[j+len_alpha+2] += dLk_dbeta[j]/(L_YkNk+1e-323);
      }
      score[len_alpha+2+len_beta] += (1-rand_indic[k])*dLk_dsigma/(L_YkNk+1e-323);
      score[len_alpha+2+len_beta+1] += rand_indic[k]*dLk_dsigma/(L_YkNk+1e-323);
      for(j = 0; j < (len_epsW); j++){
        score[j+len_alpha+2+len_beta+2] += dLk_depsW[j]/(L_Wk+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_alpha+2+len_beta+2+len_epsW] += dLk_depsX[j]/(L_Xk+1e-323);
      }


    } else{  //main sample
      A1k0 = 0;
      A1k1 = 0;
      B1k0 = 0;
      B1k1 = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        A1k0 += GH_w[q]*exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]);
        A1k1 += GH_w[q]*exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]);
        B1k0 += GH_w[q]*exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q])*GH_z[q];
        B1k1 += GH_w[q]*exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q])*GH_z[q];
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          A2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          A2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          B2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]))*GH_z[q];
          B2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]))*GH_z[q];
        }
      }

      // at this point we have A1 B1 A2 B2, eta, zeta

      // derivative chunks for size model
      if (condSize == FALSE) {
        dOmegak_dgamma0 = -sigma0[k]*B1k0/(A1k0+1e-323);
        dOmegak_dgamma1 = -sigma1[k]*B1k1/(A1k1+1e-323);
        dOmegak_dsigma0 = -gamma0[k]*B1k0/(A1k0+1e-323);
        dOmegak_dsigma1 = -gamma1[k]*B1k1/(A1k1+1e-323);
        dOmegak_dzetak0 = exp(zeta0[k])/(A1k0+1e-323);
        dOmegak_dzetak1 = exp(zeta1[k])/(A1k1+1e-323);
      } else{  // much simpler if conditional model is used
        dOmegak_dgamma0 = 0;
        dOmegak_dgamma1 = 0;
        dOmegak_dsigma0 = 0;
        dOmegak_dsigma1 = 0;
        dOmegak_dzetak0 = 1;
        dOmegak_dzetak1 = 1;
      }
      for(j = (0); j < (len_alpha); j++){
        dOmegak_dalpha0[j] =  dOmegak_dzetak0*Zk0(k,j);
        dOmegak_dalpha1[j] =  dOmegak_dzetak1*Zk1(k,j);
      }

      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = expit(eta0[i])*(1-expit(eta0[i]))/(A2ki0[i]+1e-323);
          ddeltaki_detaki1[i] = expit(eta1[i])*(1-expit(eta1[i]))/(A2ki1[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = -B2ki0[i]/(A2ki0[i]+1e-323);
          ddeltaki_sigma1[i] = -B2ki1[i]/(A2ki1[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = 1;
          ddeltaki_detaki1[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = 0;
          ddeltaki_sigma1[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_YkNk0 = 0;
      L_YkNk1 = 0;
      dLk_dalpha0 = 0*dLk_dalpha0;
      dLk_dalpha1 = 0*dLk_dalpha1;
      dLk_dgamma0 = 0;
      dLk_dgamma1 = 0;
      dLk_dbeta0 = 0*dLk_dbeta0;
      dLk_dbeta1 = 0*dLk_dbeta1;
      dLk_dsigma0 = 0;
      dLk_dsigma1 = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        // likelihood and derivative contributions form size
        L_Nk0 = exp(Nk_out[k]*(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q])-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
        L_Nk1 = exp(Nk_out[k]*(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q])-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood

        for(j = (0); j < (len_alpha); j++){
          dLNk_dalpha0[j] =  L_Nk0*dOmegak_dalpha0[j]*(Nk_out[k]-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]));
          dLNk_dalpha1[j] =  L_Nk1*dOmegak_dalpha1[j]*(Nk_out[k]-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]));
        }
        dLNk_dgamma0 = L_Nk0*(dOmegak_dgamma0+sigma0[k]*GH_z[q])*(Nk_out[k]-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]));
        dLNk_dgamma1 = L_Nk1*(dOmegak_dgamma1+sigma1[k]*GH_z[q])*(Nk_out[k]-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]));
        dLNk_dsigma0 = L_Nk0*(dOmegak_dsigma0+gamma0[k]*GH_z[q])*(Nk_out[k]-exp(Omega0[k]+gamma0[k]*sigma0[k]*GH_z[q]));
        dLNk_dsigma1 = L_Nk1*(dOmegak_dsigma1+gamma1[k]*GH_z[q])*(Nk_out[k]-exp(Omega1[k]+gamma1[k]*sigma1[k]*GH_z[q]));

        // likelihood and derivative contributions from outcome
        L_Yk0 = 1; // contribution of Yk to likelihood
        L_Yk1 = 1;
        dLYk_dsigma0 = 0;
        dLYk_dsigma1 = 0;
        dLYk_dbeta0 = 0*dLYk_dbeta0;
        dLYk_dbeta1 = 0*dLYk_dbeta1;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          L_Yk0 *=  exp(Yki[i]*(Delta0[i]+sigma0[k]*GH_z[q]))/(1+exp(Delta0[i]+sigma0[k]*GH_z[q]));
          L_Yk1 *=  exp(Yki[i]*(Delta1[i]+sigma1[k]*GH_z[q]))/(1+exp(Delta1[i]+sigma1[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta0[j] += ddeltaki_dbeta0(i,j)*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
            dLYk_dbeta1[j] += ddeltaki_dbeta1(i,j)*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          }
          dLYk_dsigma0 += (ddeltaki_sigma0[i]+GH_z[q])*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          dLYk_dsigma1 += (ddeltaki_sigma1[i]+GH_z[q])*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta0[j] =  L_Yk0*dLYk_dbeta0[j];  // scaling it by L_Yk as required
          dLYk_dbeta1[j] =  L_Yk1*dLYk_dbeta1[j];
        }
        dLYk_dsigma0 = L_Yk0*dLYk_dsigma0; // scaling it by L_Yk as required
        dLYk_dsigma1 = L_Yk1*dLYk_dsigma1;

        L_YkNk0 += GH_w[q]*L_Nk0*L_Yk0;
        L_YkNk1 += GH_w[q]*L_Nk1*L_Yk1;
        // compute main derivs here
        for(j = (0); j < (len_alpha); j++){
          dLk_dalpha0[j] += GH_w[q]*L_Yk0*dLNk_dalpha0[j];
          dLk_dalpha1[j] += GH_w[q]*L_Yk1*dLNk_dalpha1[j];
        }
        dLk_dgamma0 += GH_w[q]*L_Yk0*dLNk_dgamma0;
        dLk_dgamma1 += GH_w[q]*L_Yk1*dLNk_dgamma1;
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta0[j] += GH_w[q]*dLYk_dbeta0[j]*L_Nk0;
          dLk_dbeta1[j] += GH_w[q]*dLYk_dbeta1[j]*L_Nk1;
        }
        dLk_dsigma0 += GH_w[q]*sigma0[k]*(dLYk_dsigma0*L_Nk0+L_Yk0*dLNk_dsigma0);
        dLk_dsigma1 += GH_w[q]*sigma1[k]*(dLYk_dsigma1*L_Nk1+L_Yk1*dLNk_dsigma1);

      }
      L_Wk0 = exp(Wk[k]*(etaW0[k]))/(1+exp(etaW0[k])); // misclassification likelihood;
      L_Wk1 = exp(Wk[k]*(etaW1[k]))/(1+exp(etaW1[k]));
      L_Xk0 = exp(0*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      L_Xk1 = exp(1*(etaX[k]))/(1+exp(etaX[k]));
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW0[j] = L_Wk0*Zwk0(k,j)*(Wk[k]-expit(etaW0[k]));
        dLk_depsW1[j] = L_Wk1*Zwk1(k,j)*(Wk[k]-expit(etaW1[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX0[j] = L_Xk0*Zxk(k,j)*(0-expit(etaX[k]));
        dLk_depsX1[j] = L_Xk1*Zxk(k,j)*(1-expit(etaX[k]));
      }

      // compute denominator
      L_k = L_Wk0*L_YkNk0*L_Xk0+L_Wk1*L_YkNk1*L_Xk1;

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
        score[j] += (L_Wk0*dLk_dalpha0[j]*L_Xk0+L_Wk1*dLk_dalpha1[j]*L_Xk1)/(L_k+1e-323);
      }
      score[len_alpha] += (L_Wk0*dLk_dgamma0*L_Xk0+L_Wk1*0*L_Xk1)/(L_k+1e-323);
      score[len_alpha+1] += (L_Wk0*0*L_Xk0+L_Wk1*dLk_dgamma1*L_Xk1)/(L_k+1e-323);
      for(j = 0; j < (len_beta); j++){
        score[j+len_alpha+2] += (L_Wk0*dLk_dbeta0[j]*L_Xk0+L_Wk1*dLk_dbeta1[j]*L_Xk1)/(L_k+1e-323);
      }
      score[len_alpha+2+len_beta] += (L_Wk0*dLk_dsigma0*L_Xk0+L_Wk1*0*L_Xk1)/(L_k+1e-323);
      score[len_alpha+2+len_beta+1] += (L_Wk0*0*L_Xk0+L_Wk1*dLk_dsigma1*L_Xk1)/(L_k+1e-323);
      for(j = 0; j < (len_epsW); j++){
        score[j+len_alpha+2+len_beta+2] += (dLk_depsW0[j]*L_YkNk0*L_Xk0+dLk_depsW1[j]*L_YkNk1*L_Xk1)/(L_k+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_alpha+2+len_beta+2+len_epsW] += (L_Wk0*L_YkNk0*dLk_depsX0[j]+L_Wk1*L_YkNk1*dLk_depsX1[j])/(L_k+1e-323);
      }

    }


    cumNk += Nk_actual[k];
  }

  return score;
}





// // // // // // // // // // // // // // // // // // // // // // // //
// GLMM only


// [[Rcpp::export]]
double compute_logLik_MISclassGLMM(NumericVector sigma,
                                   NumericVector sigma0,
                                   NumericVector sigma1,
                                   NumericVector Delta,
                                   NumericVector Delta0,
                                   NumericVector Delta1,
                                   NumericVector etaW,
                                   NumericVector etaW0,
                                   NumericVector etaW1,
                                   NumericVector etaX,
                                   NumericVector Nk_actual,
                                   NumericVector Yki,
                                   NumericVector Wk,
                                   NumericVector Xk,
                                   NumericVector IDk,
                                   LogicalVector VALIDk,
                                   NumericVector GH_w,
                                   NumericVector GH_z) {

  int maxK = max(IDk); // number of clusters
  int nquad = GH_w.size(); // number of quadrature points
  int cumNk = 0; // to only subset outcomes from cluster k
  long double logL = 0; // log Likelihood (objective function to be maximized)
  int k, q, i;
  long double L_k, L_yk, L_Yk, L_yk0, L_Yk0, L_yk1, L_Yk1;
  long double L_Wk, L_Wk0, L_Wk1;
  long double L_Xk, P_Xk1;


  for(k = 0; k < maxK; k++) { // loop over clusters
    L_yk = 0; // kth contribution to the (non-log) likelihood
    L_yk0 = 0;
    L_yk1 = 0;

    for(q = 0; q < nquad; q++){ // loop over GH zeros
      // contribution of Yk to likelihood
      L_Yk = 1;
      L_Yk0 = 1;
      L_Yk1 = 1;
      for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
        L_Yk *=  exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q]));
        L_Yk0 *=  exp(Yki[i]*(Delta0[i]+sigma0[k]*GH_z[q]))/(1+exp(Delta0[i]+sigma0[k]*GH_z[q]));
        L_Yk1 *=  exp(Yki[i]*(Delta1[i]+sigma1[k]*GH_z[q]))/(1+exp(Delta1[i]+sigma1[k]*GH_z[q]));
      }
      L_yk +=  exp(log(GH_w[q])+log(L_Yk));
      L_yk0 +=  exp(log(GH_w[q])+log(L_Yk0));
      L_yk1 +=  exp(log(GH_w[q])+log(L_Yk1));
    }
    L_Wk =  exp(Wk[k]*(etaW[k]))/(1+exp(etaW[k])); // misclassification likelihood
    L_Wk0 =  exp(Wk[k]*(etaW0[k]))/(1+exp(etaW0[k]));
    L_Wk1 =  exp(Wk[k]*(etaW1[k]))/(1+exp(etaW1[k]));
    // single exposure likelihood
    L_Xk =  exp(Xk[k]*(etaX[k]))/(1+exp(etaX[k]));
    P_Xk1 =  exp(1*(etaX[k]))/(1+exp(etaX[k]));

    if(VALIDk[k] == TRUE){
      L_k = L_Wk*L_yk*L_Xk;
    } else {
      L_k = L_Wk0*L_yk0*(1-P_Xk1)+L_Wk1*L_yk1*P_Xk1;
    }

    cumNk += Nk_actual[k];

    logL += log(L_k+1e-323);
  }

  return logL;
}



// note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// [[Rcpp::export]]
NumericVector compute_score_MISclassGLMM(NumericVector sigma,
                                         NumericVector sigma0,
                                         NumericVector sigma1,
                                         NumericVector Delta,
                                         NumericVector Delta0,
                                         NumericVector Delta1,
                                         NumericVector eta,
                                         NumericVector eta0,
                                         NumericVector eta1,
                                         NumericVector etaW,
                                         NumericVector etaW0,
                                         NumericVector etaW1,
                                         NumericVector etaX,
                                         NumericVector Nk_actual,
                                         NumericVector Yki,
                                         NumericMatrix Xki,
                                         NumericMatrix Xki0,
                                         NumericMatrix Xki1,
                                         NumericVector Wk,
                                         NumericMatrix Zwk,
                                         NumericMatrix Zwk0,
                                         NumericMatrix Zwk1,
                                         NumericVector Xk,
                                         NumericMatrix Zxk,
                                         NumericVector IDk,
                                         LogicalVector VALIDk,
                                         bool condOut,
                                         NumericVector GH_w,
                                         NumericVector GH_z,
                                         int len_beta,
                                         int len_epsW,
                                         int len_epsX) {


  int maxK = Nk_actual.size(); // number of clusters
  int N = Delta.size(); // number of obs
  int nquad = GH_w.size(); // number of quadrature points
  int cumNk = 0; // to only subset outcomes from cluster k
  int k, q, i, j;
  NumericVector score(len_beta+1+len_epsW+len_epsX);
  long double L_k, L_yk, L_yk0, L_yk1, L_Yk, L_k0, L_Yk0, L_k1, L_Yk1;
  long double L_Wk, L_Wk0, L_Wk1, L_Xk, L_Xk0, L_Xk1;
  NumericVector A2ki(N), B2ki(N), A2ki0(N), B2ki0(N), A2ki1(N), B2ki1(N);
  NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N), ddeltaki_detaki0(N), ddeltaki_sigma0(N), ddeltaki_detaki1(N), ddeltaki_sigma1(N);
  NumericMatrix ddeltaki_dbeta(N,len_beta), ddeltaki_dbeta0(N,len_beta), ddeltaki_dbeta1(N,len_beta);
  long double dLYk_dsigma, dLk_dsigma, dLYk_dsigma0, dLk_dsigma0, dLYk_dsigma1, dLk_dsigma1;
  NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta), dLYk_dbeta0(len_beta), dLk_dbeta0(len_beta), dLYk_dbeta1(len_beta), dLk_dbeta1(len_beta);
  NumericVector dLk_depsW(len_epsW), dLk_depsX(len_epsX), dLk_depsW0(len_epsW), dLk_depsX0(len_epsX), dLk_depsW1(len_epsW), dLk_depsX1(len_epsX);

  for(k = 0; k < maxK; k++) { // loop over clusters

    if(VALIDk[k] == TRUE){
      for(q = 0; q < nquad; q++){ // loop over GH zeros

        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes

          A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));

          B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];

        }

      }

      // at this point we have A2 B2, eta

      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_yk = 0;
      dLk_dbeta = 0*dLk_dbeta;
      dLk_dsigma = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros

        // likelihood and derivative contributions from outcome
        L_Yk = 1; // contribution of Yk to likelihood
        dLYk_dsigma = 0;
        dLYk_dbeta = 0*dLYk_dbeta;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes

          L_Yk *=  exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta[j] += ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
          }
          dLYk_dsigma += (ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
        }
        dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required

        L_yk += GH_w[q]*L_Yk;
        // compute main derivs here
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j];
        }
        dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma);
      }

      L_Wk = exp(Wk[k]*(etaW[k]))/(1+exp(etaW[k])); // misclassification likelihood;
      L_Xk = exp(Xk[k]*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW[j] = L_Wk*Zwk(k,j)*(Wk[k]-expit(etaW[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX[j] = L_Xk*Zxk(k,j)*(Xk[k]-expit(etaX[k]));
      }

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_beta); j++){
        score[j] += dLk_dbeta[j]/(L_yk+1e-323);
      }
      score[len_beta] += dLk_dsigma/(L_yk+1e-323);
      for(j = 0; j < (len_epsW); j++){
        score[j+len_beta+1] += dLk_depsW[j]/(L_Wk+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_beta+1+len_epsW] += dLk_depsX[j]/(L_Xk+1e-323);
      }
    } else{ //main sample

      for(q = 0; q < nquad; q++){ // loop over GH zeros
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          A2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          A2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          B2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]))*GH_z[q];
          B2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]))*GH_z[q];
        }
      }

      // at this point we have A2 B2, eta
      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = expit(eta0[i])*(1-expit(eta0[i]))/(A2ki0[i]+1e-323);
          ddeltaki_detaki1[i] = expit(eta1[i])*(1-expit(eta1[i]))/(A2ki1[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = -B2ki0[i]/(A2ki0[i]+1e-323);
          ddeltaki_sigma1[i] = -B2ki1[i]/(A2ki1[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = 1;
          ddeltaki_detaki1[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = 0;
          ddeltaki_sigma1[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_yk0 = 0;
      L_yk1 = 0;
      dLk_dbeta0 = 0*dLk_dbeta0;
      dLk_dbeta1 = 0*dLk_dbeta1;
      dLk_dsigma0 = 0;
      dLk_dsigma1 = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        // likelihood and derivative contributions from outcome
        L_Yk0 = 1; // contribution of Yk to likelihood
        L_Yk1 = 1;
        dLYk_dsigma0 = 0;
        dLYk_dsigma1 = 0;
        dLYk_dbeta0 = 0*dLYk_dbeta0;
        dLYk_dbeta1 = 0*dLYk_dbeta1;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          L_Yk0 *=  exp(Yki[i]*(Delta0[i]+sigma0[k]*GH_z[q]))/(1+exp(Delta0[i]+sigma0[k]*GH_z[q]));
          L_Yk1 *=  exp(Yki[i]*(Delta1[i]+sigma1[k]*GH_z[q]))/(1+exp(Delta1[i]+sigma1[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta0[j] += ddeltaki_dbeta0(i,j)*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
            dLYk_dbeta1[j] += ddeltaki_dbeta1(i,j)*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          }
          dLYk_dsigma0 += (ddeltaki_sigma0[i]+GH_z[q])*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          dLYk_dsigma1 += (ddeltaki_sigma1[i]+GH_z[q])*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta0[j] =  L_Yk0*dLYk_dbeta0[j];  // scaling it by L_Yk as required
          dLYk_dbeta1[j] =  L_Yk1*dLYk_dbeta1[j];
        }
        dLYk_dsigma0 = L_Yk0*dLYk_dsigma0; // scaling it by L_Yk as required
        dLYk_dsigma1 = L_Yk1*dLYk_dsigma1;

        L_yk0 += GH_w[q]*L_Yk0;
        L_yk1 += GH_w[q]*L_Yk1;
        // compute main derivs here
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta0[j] += GH_w[q]*dLYk_dbeta0[j];
          dLk_dbeta1[j] += GH_w[q]*dLYk_dbeta1[j];
        }
        dLk_dsigma0 += GH_w[q]*sigma0[k]*(dLYk_dsigma0);
        dLk_dsigma1 += GH_w[q]*sigma1[k]*(dLYk_dsigma1);
      }

      L_Wk0 = exp(Wk[k]*(etaW0[k]))/(1+exp(etaW0[k])); // misclassification likelihood;
      L_Wk1 = exp(Wk[k]*(etaW1[k]))/(1+exp(etaW1[k]));
      L_Xk0 = exp(0*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      L_Xk1 = exp(1*(etaX[k]))/(1+exp(etaX[k]));
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW0[j] = L_Wk0*Zwk0(k,j)*(Wk[k]-expit(etaW0[k]));
        dLk_depsW1[j] = L_Wk1*Zwk1(k,j)*(Wk[k]-expit(etaW1[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX0[j] = L_Xk0*Zxk(k,j)*(0-expit(etaX[k]));
        dLk_depsX1[j] = L_Xk1*Zxk(k,j)*(1-expit(etaX[k]));
      }

      // compute denominator
      L_k = L_Wk0*L_yk0*L_Xk0+L_Wk1*L_yk1*L_Xk1;

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_beta); j++){
        score[j] += (L_Wk0*dLk_dbeta0[j]*L_Xk0+L_Wk1*dLk_dbeta1[j]*L_Xk1)/(L_k+1e-323);
      }
      score[len_beta] += (L_Wk0*dLk_dsigma0*L_Xk0+L_Wk1*dLk_dsigma1*L_Xk1)/(L_k+1e-323);
      for(j = 0; j < (len_epsW); j++){
        score[j+len_beta+1] += (dLk_depsW0[j]*L_yk0*L_Xk0+dLk_depsW1[j]*L_yk1*L_Xk1)/(L_k+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_beta+1+len_epsW] += (L_Wk0*L_yk0*dLk_depsX0[j]+L_Wk1*L_yk1*dLk_depsX1[j])/(L_k+1e-323);
      }

    }

    cumNk += Nk_actual[k];
  }

  return score;
}


// note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// [[Rcpp::export]]
NumericVector compute_score_MISclassGLMM_slopes(NumericVector sigma,
                                                NumericVector sigma0,
                                                NumericVector sigma1,
                                                NumericVector Delta,
                                                NumericVector Delta0,
                                                NumericVector Delta1,
                                                NumericVector eta,
                                                NumericVector eta0,
                                                NumericVector eta1,
                                                NumericVector etaW,
                                                NumericVector etaW0,
                                                NumericVector etaW1,
                                                NumericVector etaX,
                                                NumericVector rand_indic,
                                                NumericVector Nk_actual,
                                                NumericVector Yki,
                                                NumericMatrix Xki,
                                                NumericMatrix Xki0,
                                                NumericMatrix Xki1,
                                                NumericVector Wk,
                                                NumericMatrix Zwk,
                                                NumericMatrix Zwk0,
                                                NumericMatrix Zwk1,
                                                NumericVector Xk,
                                                NumericMatrix Zxk,
                                                NumericVector IDk,
                                                LogicalVector VALIDk,
                                                bool condOut,
                                                NumericVector GH_w,
                                                NumericVector GH_z,
                                                int len_beta,
                                                int len_epsW,
                                                int len_epsX) {


  int maxK = Nk_actual.size(); // number of clusters
  int N = Delta.size(); // number of obs
  int nquad = GH_w.size(); // number of quadrature points
  int cumNk = 0; // to only subset outcomes from cluster k
  int k, q, i, j;
  NumericVector score(len_beta+2+len_epsW+len_epsX); //2 sigmas
  long double L_k;
  long double L_yk, L_Yk, L_yk0, L_Yk0, L_yk1, L_Yk1;
  long double L_Wk, L_Xk, L_Wk0, L_Xk0, L_Wk1, L_Xk1;
  NumericVector A2ki(N), B2ki(N), A2ki0(N), B2ki0(N), A2ki1(N), B2ki1(N);
  NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N), ddeltaki_detaki0(N), ddeltaki_sigma0(N), ddeltaki_detaki1(N), ddeltaki_sigma1(N);
  NumericMatrix ddeltaki_dbeta(N,len_beta), ddeltaki_dbeta0(N,len_beta), ddeltaki_dbeta1(N,len_beta);
  long double dLYk_dsigma, dLk_dsigma, dLYk_dsigma0, dLk_dsigma0, dLYk_dsigma1, dLk_dsigma1;
  NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta), dLYk_dbeta0(len_beta), dLk_dbeta0(len_beta), dLYk_dbeta1(len_beta), dLk_dbeta1(len_beta);
  NumericVector dLk_depsW(len_epsW), dLk_depsX(len_epsX), dLk_depsW0(len_epsW), dLk_depsX0(len_epsX), dLk_depsW1(len_epsW), dLk_depsX1(len_epsX);

  for(k = 0; k < maxK; k++) { // loop over clusters

    if(VALIDk[k] == TRUE){

      for(q = 0; q < nquad; q++){ // loop over GH zeros

        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes

          A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));

          B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];

        }

      }

      // at this point we have A2 B2, eta

      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
          }
          ddeltaki_sigma[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_yk = 0;
      dLk_dbeta = 0*dLk_dbeta;
      dLk_dsigma = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros

        // likelihood and derivative contributions from outcome
        L_Yk = 1; // contribution of Yk to likelihood
        dLYk_dsigma = 0;
        dLYk_dbeta = 0*dLYk_dbeta;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes

          L_Yk *=  exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta[j] += ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
          }
          dLYk_dsigma += (ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
        }
        dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required



        L_yk += GH_w[q]*L_Yk;
        // compute main derivs here
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j];
        }
        dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma);


      }

      L_Wk = exp(Wk[k]*(etaW[k]))/(1+exp(etaW[k])); // misclassification likelihood;
      L_Xk = exp(Xk[k]*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW[j] = L_Wk*Zwk(k,j)*(Wk[k]-expit(etaW[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX[j] = L_Xk*Zxk(k,j)*(Xk[k]-expit(etaX[k]));
      }

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_beta); j++){
        score[j] += dLk_dbeta[j]/(L_yk+1e-323);
      }
      score[len_beta] += (1-rand_indic[k])*dLk_dsigma/(L_yk+1e-323);
      score[len_beta+1] += rand_indic[k]*dLk_dsigma/(L_yk+1e-323);

      for(j = 0; j < (len_epsW); j++){
        score[j+len_beta+2] += dLk_depsW[j]/(L_Wk+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_beta+2+len_epsW] += dLk_depsX[j]/(L_Xk+1e-323);
      }



    } else{ // main sample

      for(q = 0; q < nquad; q++){ // loop over GH zeros
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          A2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          A2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          B2ki0[i] += GH_w[q]*expit(Delta0[i]+sigma0[k]*GH_z[q])*(1-expit(Delta0[i]+sigma0[k]*GH_z[q]))*GH_z[q];
          B2ki1[i] += GH_w[q]*expit(Delta1[i]+sigma1[k]*GH_z[q])*(1-expit(Delta1[i]+sigma1[k]*GH_z[q]))*GH_z[q];
        }
      }

      // at this point we have A2 B2, eta
      // derivative chunks for outcome model
      if (condOut == FALSE) {
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = expit(eta0[i])*(1-expit(eta0[i]))/(A2ki0[i]+1e-323);
          ddeltaki_detaki1[i] = expit(eta1[i])*(1-expit(eta1[i]))/(A2ki1[i]+1e-323);
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = -B2ki0[i]/(A2ki0[i]+1e-323);
          ddeltaki_sigma1[i] = -B2ki1[i]/(A2ki1[i]+1e-323);
        }
      } else { // much simpler if conditional model is used
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          ddeltaki_detaki0[i] = 1;
          ddeltaki_detaki1[i] = 1;
          for(j = (0); j < (len_beta); j++){
            ddeltaki_dbeta0(i,j) =  ddeltaki_detaki0[i]*Xki0(i,j);
            ddeltaki_dbeta1(i,j) =  ddeltaki_detaki1[i]*Xki1(i,j);
          }
          ddeltaki_sigma0[i] = 0;
          ddeltaki_sigma1[i] = 0;
        }
      }

      // at this point we have deriv chunks that will be needed for the next round of integration
      L_yk0 = 0;
      L_yk1 = 0;
      dLk_dbeta0 = 0*dLk_dbeta0;
      dLk_dbeta1 = 0*dLk_dbeta1;
      dLk_dsigma0 = 0;
      dLk_dsigma1 = 0;
      for(q = 0; q < nquad; q++){ // loop over GH zeros
        // likelihood and derivative contributions from outcome
        L_Yk0 = 1; // contribution of Yk to likelihood
        L_Yk1 = 1;
        dLYk_dsigma0 = 0;
        dLYk_dsigma1 = 0;
        dLYk_dbeta0 = 0*dLYk_dbeta0;
        dLYk_dbeta1 = 0*dLYk_dbeta1;
        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
          L_Yk0 *=  exp(Yki[i]*(Delta0[i]+sigma0[k]*GH_z[q]))/(1+exp(Delta0[i]+sigma0[k]*GH_z[q]));
          L_Yk1 *=  exp(Yki[i]*(Delta1[i]+sigma1[k]*GH_z[q]))/(1+exp(Delta1[i]+sigma1[k]*GH_z[q]));
          for(j = (0); j < (len_beta); j++){
            dLYk_dbeta0[j] += ddeltaki_dbeta0(i,j)*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
            dLYk_dbeta1[j] += ddeltaki_dbeta1(i,j)*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
          }
          dLYk_dsigma0 += (ddeltaki_sigma0[i]+GH_z[q])*(Yki[i]-expit(Delta0[i]+sigma0[k]*GH_z[q]));
          dLYk_dsigma1 += (ddeltaki_sigma1[i]+GH_z[q])*(Yki[i]-expit(Delta1[i]+sigma1[k]*GH_z[q]));
        }
        for(j = (0); j < (len_beta); j++){
          dLYk_dbeta0[j] =  L_Yk0*dLYk_dbeta0[j];  // scaling it by L_Yk as required
          dLYk_dbeta1[j] =  L_Yk1*dLYk_dbeta1[j];
        }
        dLYk_dsigma0 = L_Yk0*dLYk_dsigma0; // scaling it by L_Yk as required
        dLYk_dsigma1 = L_Yk1*dLYk_dsigma1;


        L_yk0 += GH_w[q]*L_Yk0;
        L_yk1 += GH_w[q]*L_Yk1;
        // compute main derivs here
        for(j = (0); j < (len_beta); j++){
          dLk_dbeta0[j] += GH_w[q]*dLYk_dbeta0[j];
          dLk_dbeta1[j] += GH_w[q]*dLYk_dbeta1[j];
        }
        dLk_dsigma0 += GH_w[q]*sigma0[k]*(dLYk_dsigma0);
        dLk_dsigma1 += GH_w[q]*sigma1[k]*(dLYk_dsigma1);
      }

      L_Wk0 = exp(Wk[k]*(etaW0[k]))/(1+exp(etaW0[k])); // misclassification likelihood;
      L_Wk1 = exp(Wk[k]*(etaW1[k]))/(1+exp(etaW1[k]));
      L_Xk0 = exp(0*(etaX[k]))/(1+exp(etaX[k])); //exposure likelihood
      L_Xk1 = exp(1*(etaX[k]))/(1+exp(etaX[k]));
      for(j = 0; j < (len_epsW); j++){
        dLk_depsW0[j] = L_Wk0*Zwk0(k,j)*(Wk[k]-expit(etaW0[k]));
        dLk_depsW1[j] = L_Wk1*Zwk1(k,j)*(Wk[k]-expit(etaW1[k]));
      }
      for(j = 0; j < (len_epsX); j++){
        dLk_depsX0[j] = L_Xk0*Zxk(k,j)*(0-expit(etaX[k]));
        dLk_depsX1[j] = L_Xk1*Zxk(k,j)*(1-expit(etaX[k]));
      }

      // compute denominator
      L_k = L_Wk0*L_yk0*L_Xk0+L_Wk1*L_yk1*L_Xk1;

      //then scale all by L_k and sum over k
      for(j = 0; j < (len_beta); j++){
        score[j] += (L_Wk0*dLk_dbeta0[j]*L_Xk0+L_Wk1*dLk_dbeta1[j]*L_Xk1)/(L_k+1e-323);
      }
      score[len_beta] += (L_Wk0*dLk_dsigma0*L_Xk0)/(L_k+1e-323);
      score[len_beta+1] += (L_Wk1*dLk_dsigma1*L_Xk1)/(L_k+1e-323);
      for(j = 0; j < (len_epsW); j++){
        score[j+len_beta+2] += (dLk_depsW0[j]*L_yk0*L_Xk0+dLk_depsW1[j]*L_yk1*L_Xk1)/(L_k+1e-323);
      }
      for(j = 0; j < (len_epsX); j++){
        score[j+len_beta+2+len_epsW] += (L_Wk0*L_yk0*dLk_depsX0[j]+L_Wk1*L_yk1*dLk_depsX1[j])/(L_k+1e-323);
      }

    }

    cumNk += Nk_actual[k];
  }

  return score;
}



// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_JMMICS(NumericVector gamma,
//                                     NumericVector sigma,
//                                    NumericVector Omega,
//                                    NumericVector Delta,
//                                    NumericVector zeta,
//                                    NumericVector eta,
//                                    NumericVector Nk,
//                                    NumericVector Nk_actual,
//                                    NumericMatrix Zk,
//                                    NumericVector Yki,
//                                    NumericMatrix Xki,
//                                    NumericVector IDk,
//                                    NumericVector Wki,
//                                    int minNk,
//                                    bool condSize,
//                                    bool condOut,
//                                    NumericVector GH_w,
//                                    NumericVector GH_z,
//                                    int len_alpha,
//                                    int len_beta) {
//
//         int maxK = Omega.size(); // number of clusters
//         int N = Delta.size(); // number of obs
//         int nquad = GH_w.size(); // number of quadrature points
//         int cumNk = 0; // to only subset outcomes from cluster k
//         int k, q, i, j;
//         NumericVector Nk_out = Nk - minNk;
//         NumericVector logNk_factorial = lgamma(Nk_out+1);
//         NumericVector score(len_alpha+len_beta+1+1);
//         long double zetak, L_k, L_Nk, L_Yk, A1k, B1k;
//         NumericVector A2ki(N), B2ki(N);
//         long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//         NumericVector dOmegak_dalpha(len_alpha);
//         NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//         NumericMatrix ddeltaki_dbeta(N,len_beta);
//         long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma;
//         NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha);
//         NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//         NumericMatrix cheese(len_alpha+len_beta+1+1,len_alpha+len_beta+1+1);
//
//         for(k = 0; k < maxK; k++) { // loop over clusters
//
//           zetak = zeta[k];
//
//           A1k = 0;
//           B1k = 0;
//           for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//             A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//             B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//             for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//               A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//               B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//             }
//
//           }
//
//           // at this point we have A1 B1 A2 B2, eta, zeta
//
//           // derivative chunks for size model
//           if (condSize == FALSE) {
//             dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//             dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//             dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//           } else{  // much simpler if conditional model is used
//             dOmegak_dgamma = 0;
//             dOmegak_dsigma = 0;
//             dOmegak_dzetak = 1;
//           }
//           for(j = (0); j < (len_alpha); j++){
//             dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//           }
//
//           // derivative chunks for outcome model
//           if (condOut == FALSE) {
//             for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//               ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//               for(j = (0); j < (len_beta); j++){
//                 ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//               }
//               ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//             }
//           } else { // much simpler if conditional model is used
//             for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//               ddeltaki_detaki[i] = 1;
//                 for(j = (0); j < (len_beta); j++){
//                   ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//                 }
//               ddeltaki_sigma[i] = 0;
//             }
//           }
//
//           // at this point we have deriv chunks that will be needed for the next round of integration
//           L_k = 0;
//           dLk_dalpha = 0*dLk_dalpha;
//           dLk_dgamma = 0;
//           dLk_dbeta = 0*dLk_dbeta;
//           dLk_dsigma = 0;
//           for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//             // likelihood and derivative contributions form size
//             L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//             for(j = (0); j < (len_alpha); j++){
//               dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//             }
//             dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//             dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//
//
//           // likelihood and derivative contributions from outcome
//             L_Yk = 1; // contribution of Yk to likelihood
//             dLYk_dsigma = 0;
//             dLYk_dbeta = 0*dLYk_dbeta;
//             for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//               L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//               for(j = (0); j < (len_beta); j++){
//                 dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//               }
//               dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//             }
//             for(j = (0); j < (len_beta); j++){
//               dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//             }
//             dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//             L_k += GH_w[q]*L_Nk*L_Yk;
//             // compute main derivs here
//             for(j = (0); j < (len_alpha); j++){
//               dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
//             }
//             dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
//             for(j = (0); j < (len_beta); j++){
//               dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
//             }
//             dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma);
//
//
//
//           }
//
//           //then scale all by L_k and sum over k
//           for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//             score[j] = dLk_dalpha[j]/(L_k+1e-323);
//           }
//           score[len_alpha] = dLk_dgamma/(L_k+1e-323);
//           for(j = 0; j < (len_beta); j++){
//             score[j+len_alpha+1] = dLk_dbeta[j]/(L_k+1e-323);
//           }
//           score[len_alpha+1+len_beta] = dLk_dsigma/(L_k+1e-323);
//
//
//           cheese += square_vect(score);
//
//           cumNk += Nk_actual[k];
//         }
//
//         return cheese;
// }
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_JMMICS_slopes(NumericVector gamma,
//                                            NumericVector sigma,
//                                            NumericVector Omega,
//                                            NumericVector Delta,
//                                            NumericVector zeta,
//                                            NumericVector eta,
//                                            NumericVector rand_indic,
//                                              NumericVector Nk,
//                                              NumericVector Nk_actual,
//                                              NumericMatrix Zk,
//                                              NumericVector Yki,
//                                              NumericMatrix Xki,
//                                              NumericVector IDk,
//                                              NumericVector Wki,
//                                              int minNk,
//                                              bool condSize,
//                                              bool condOut,
//                                              NumericVector GH_w,
//                                              NumericVector GH_z,
//                                              int len_alpha,
//                                              int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector score(len_alpha+len_beta+2+2); // 2 gammas and 2 sigmas
//   long double zetak, L_k, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma;
//   NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//   NumericMatrix cheese(len_alpha+len_beta+2+2,len_alpha+len_beta+2+2); // 2 gammas and 2 sigmas
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     zetak = zeta[k];
//
//     A1k = 0;
//     B1k = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2, eta, zeta
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//     dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dalpha = 0*dLk_dalpha;
//     dLk_dgamma = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//     // likelihood and derivative contributions form size
//     L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//     for(j = (0); j < (len_alpha); j++){
//       dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//     }
//     dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//     dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//
//
//     // likelihood and derivative contributions from outcome
//     L_Yk = 1; // contribution of Yk to likelihood
//     dLYk_dsigma = 0;
//     dLYk_dbeta = 0*dLYk_dbeta;
//     for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//       L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//     }
//     for(j = (0); j < (len_beta); j++){
//       dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//     }
//     dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//     L_k += GH_w[q]*L_Nk*L_Yk;
//     // compute main derivs here
//     for(j = (0); j < (len_alpha); j++){
//       dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
//     }
//     dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
//     for(j = (0); j < (len_beta); j++){
//       dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
//     }
//     dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma);
//
//
//
//   }
//
//   //then scale all by L_k and sum over k
//   for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//   score[j] = dLk_dalpha[j]/(L_k+1e-323);
//   }
//   score[len_alpha] = (1-rand_indic[k])*dLk_dgamma/(L_k+1e-323);
//     score[len_alpha+1] = rand_indic[k]*dLk_dgamma/(L_k+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+2] = dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_alpha+2+len_beta] = (1-rand_indic[k])*dLk_dsigma/(L_k+1e-323);
//     score[len_alpha+2+len_beta+1] = rand_indic[k]*dLk_dsigma/(L_k+1e-323);
//
//
//
//     cheese += square_vect(score);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return cheese;
// }
//
// // [[Rcpp::export]]
// NumericVector compute_EB_JMMICS(NumericVector gamma,
//                                 NumericVector sigma,
//                          NumericVector Omega,
//                          NumericVector Delta,
//                          NumericVector Nk,
//                          NumericVector Nk_actual,
//                          NumericVector Yki,
//                          NumericVector IDk,
//                          NumericVector Wki,
//                          int minNk,
//                          NumericVector GH_w,
//                          NumericVector GH_z) {
//
//     int maxK = Omega.size(); // number of clusters
//     int nquad = GH_w.size(); // number of quadrature points
//     int cumNk = 0; // to only subset outcomes from cluster k
//     int k, q, i;
//     double L_Nk, L_Yk, L_k, zL_k;
//     NumericVector pred(maxK);
//     NumericVector Nk_out = Nk - minNk; // adjust for Poisson + minimum
//     NumericVector logNk_factorial = lgamma(Nk_out+1);
//
//     for(k = 0; k < maxK; k++) { // loop over clusters
//
//         L_k = 0; // kth denominator
//         zL_k = 0; // kth numerator
//
//         for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//             L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//
//             L_Yk = 1; // contribution of Yk to likelihood
//
//             for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//                 L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//
//             }
//
//             L_k += GH_w[q]*L_Nk*L_Yk; // denominator
//             zL_k += GH_w[q]*L_Nk*L_Yk*sigma[k]*GH_z[q]; // numerator
//         }
//
//         cumNk += Nk_actual[k];
//
//         pred[k] = zL_k/(L_k+1e-323);
//     }
//
//     return pred;
// }
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// //
// // Negative Binomial Variants of logLik, score, cheese, and EB  //
// //
//
//
// // [[Rcpp::export]]
// double compute_logLik_JMMICS_NegBin(NumericVector gamma,
//                                     NumericVector sigma,
//                                    double tau, // dispersion parameter
//                                    NumericVector Omega,
//                                    NumericVector Delta,
//                                    NumericVector Nk,
//                                    NumericVector Nk_actual,
//                                    NumericVector Yki,
//                                    NumericVector IDk,
//                                    NumericVector Wki,
//                                    int minNk,
//                                    NumericVector GH_w,
//                                    NumericVector GH_z) {
//
//   int maxK = max(IDk); // number of clusters
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   long double logL = 0; // log Likelihood (objective function to be maximized)
//   int k, q, i;
//   long double L_k, L_Nk, L_Yk; // lam;
//   NumericVector Nk_out = Nk - minNk; // adjust for Poisson + minimum
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector lgam_Nktau = lgamma(Nk_out+tau);
//   double lgam_tau = lgamma(tau);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     L_k = 0; // kth contribution to the (non-log) likelihood
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//       // its really logLNk
//       //L_Nk =  Nk_out[k]*(Omega[k]+gamma[k]*sigma*GH_z[q]) -exp(Omega[k]+gamma[k]*sigma*GH_z[q])-logNk_factorial[k];   // contribution of Nk to likelihood
//       L_Nk =  lgam_Nktau[k]-logNk_factorial[k]-lgam_tau +tau*log(tau)+Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-(tau+Nk_out[k])*log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau);   // contribution of Nk to log likelihood
//       //lam =  exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//       //L_Nk =  ::Rf_dpois(Nk_out[k],lam,1);
//
//       L_Yk = 1; // contribution of Yk to likelihood
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//
//       }
//
//       L_k +=  exp(log(GH_w[q])+L_Nk+log(L_Yk));
//     }
//
//     cumNk += Nk_actual[k];
//
//     logL += log(L_k+1e-323);
//   }
//
//   return logL;
// }
//
//
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // same is true for tau!
// // [[Rcpp::export]]
// NumericVector compute_score_JMMICS_NegBin(NumericVector gamma,
//                                           NumericVector sigma,
//                                    double tau,
//                                    double DG_tau,
//                                    NumericVector Omega,
//                                    NumericVector Delta,
//                                    NumericVector zeta,
//                                    NumericVector eta,
//                                    NumericVector Nk,
//                                    NumericVector Nk_actual,
//                                    NumericMatrix Zk,
//                                    NumericVector Yki,
//                                    NumericMatrix Xki,
//                                    NumericVector IDk,
//                                    NumericVector Wki,
//                                    int minNk,
//                                    bool condSize,
//                                    bool condOut,
//                                    NumericVector GH_w,
//                                    NumericVector GH_z,
//                                    int len_alpha,
//                                    int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector lgam_Nktau = lgamma(Nk_out+tau);
//   double lgam_tau = lgamma(tau);
//   NumericVector DG_Nktau = digamma(Nk_out+tau);
//   NumericVector score(len_alpha+len_beta+1+1+1); // +1 for tau
//   long double zetak, L_k, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma, dLNk_dtau, dLk_dtau;
//   NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     zetak = zeta[k];
//
//     A1k = 0;
//     B1k = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2, eta, zeta
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//     dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dalpha = 0*dLk_dalpha;
//     dLk_dgamma = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     dLk_dtau = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//     // likelihood and derivative contributions form size
//     //L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//     L_Nk = exp(lgam_Nktau[k]-logNk_factorial[k]-lgam_tau+tau*log(tau)+Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-(tau+Nk_out[k])*log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//     for(j = (0); j < (len_alpha); j++){
//       //dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//     }
//     //dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//     dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//     //dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//     dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//     //
//     dLNk_dtau = L_Nk*(DG_Nktau[k]-DG_tau+log(tau)+1-log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau)-(tau+Nk_out[k])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau) );
//
//     // likelihood and derivative contributions from outcome
//     L_Yk = 1; // contribution of Yk to likelihood
//     dLYk_dsigma = 0;
//     dLYk_dbeta = 0*dLYk_dbeta;
//     for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//     L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//     }
//     for(j = (0); j < (len_beta); j++){
//       dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//     }
//     dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//     L_k += GH_w[q]*L_Nk*L_Yk;
//     // compute main derivs here
//     for(j = (0); j < (len_alpha); j++){
//       dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
//     }
//     dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
//     for(j = (0); j < (len_beta); j++){
//       dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
//     }
//     dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma); //multiply by sigma because max wrt logsigma
//     dLk_dtau += GH_w[q]*tau*L_Yk*dLNk_dtau; //multiply by tau because max wrt logtau
//
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//       score[j] += dLk_dalpha[j]/(L_k+1e-323);
//     }
//     score[len_alpha] += dLk_dgamma/(L_k+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+1] += dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_alpha+1+len_beta] += dLk_dsigma/(L_k+1e-323);
//     score[len_alpha+1+len_beta+1] += dLk_dtau/(L_k+1e-323);
//
//
//     cumNk += Nk_actual[k];
//   }
//
//   return score;
// }
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // same is true for tau!
// // [[Rcpp::export]]
// NumericVector compute_score_JMMICS_NegBin_slopes(NumericVector gamma,
//                                           NumericVector sigma,
//                                           double tau,
//                                           double DG_tau,
//                                           NumericVector Omega,
//                                           NumericVector Delta,
//                                           NumericVector zeta,
//                                           NumericVector eta,
//                                           NumericVector rand_indic,
//                                           NumericVector Nk,
//                                           NumericVector Nk_actual,
//                                           NumericMatrix Zk,
//                                           NumericVector Yki,
//                                           NumericMatrix Xki,
//                                           NumericVector IDk,
//                                           NumericVector Wki,
//                                           int minNk,
//                                           bool condSize,
//                                           bool condOut,
//                                           NumericVector GH_w,
//                                           NumericVector GH_z,
//                                           int len_alpha,
//                                           int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector lgam_Nktau = lgamma(Nk_out+tau);
//   double lgam_tau = lgamma(tau);
//   NumericVector DG_Nktau = digamma(Nk_out+tau);
//   NumericVector score(len_alpha+len_beta+2+2+1); // +1 for tau //2 gammas and 2 sigmas
//   long double zetak, L_k, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma, dLNk_dtau, dLk_dtau;
//   NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     zetak = zeta[k];
//
//     A1k = 0;
//     B1k = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2, eta, zeta
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//       dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//       // derivative chunks for outcome model
//       if (condOut == FALSE) {
//         for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//           ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//           for(j = (0); j < (len_beta); j++){
//             ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//           }
//           ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//         }
//       } else { // much simpler if conditional model is used
//         for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dalpha = 0*dLk_dalpha;
//     dLk_dgamma = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     dLk_dtau = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions form size
//       //L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//       L_Nk = exp(lgam_Nktau[k]-logNk_factorial[k]-lgam_tau+tau*log(tau)+Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-(tau+Nk_out[k])*log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       for(j = (0); j < (len_alpha); j++){
//         //dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//         dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       }
//       //dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       //dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       //
//       dLNk_dtau = L_Nk*(DG_Nktau[k]-DG_tau+log(tau)+1-log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau)-(tau+Nk_out[k])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau) );
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//       L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k += GH_w[q]*L_Nk*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_alpha); j++){
//         dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
//       }
//       dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
//       for(j = (0); j < (len_beta); j++){
//         dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
//       }
//       dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma); //multiply by sigma because max wrt logsigma
//       dLk_dtau += GH_w[q]*tau*L_Yk*dLNk_dtau; //multiply by tau because max wrt logtau
//
//
//
//     }
//
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//     score[j] += dLk_dalpha[j]/(L_k+1e-323);
//     }
//     score[len_alpha] += (1-rand_indic[k])*dLk_dgamma/(L_k+1e-323);
//     score[len_alpha+1] += rand_indic[k]*dLk_dgamma/(L_k+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+2] += dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_alpha+2+len_beta] += (1-rand_indic[k])*dLk_dsigma/(L_k+1e-323);
//     score[len_alpha+2+len_beta+1] += rand_indic[k]*dLk_dsigma/(L_k+1e-323);
//     score[len_alpha+2+len_beta+2] += dLk_dtau/(L_k+1e-323);
//
//
//
//
//     cumNk += Nk_actual[k];
//   }
//
//   return score;
// }
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // tau as well!
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_JMMICS_NegBin(NumericVector gamma,
//                                            NumericVector sigma,
//                                           double tau,
//                                           double DG_tau,
//                                           NumericVector Omega,
//                                           NumericVector Delta,
//                                           NumericVector zeta,
//                                           NumericVector eta,
//                                           NumericVector Nk,
//                                           NumericVector Nk_actual,
//                                           NumericMatrix Zk,
//                                           NumericVector Yki,
//                                           NumericMatrix Xki,
//                                           NumericVector IDk,
//                                           NumericVector Wki,
//                                           int minNk,
//                                           bool condSize,
//                                           bool condOut,
//                                           NumericVector GH_w,
//                                           NumericVector GH_z,
//                                           int len_alpha,
//                                           int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector lgam_Nktau = lgamma(Nk_out+tau);
//   double lgam_tau = lgamma(tau);
//   NumericVector DG_Nktau = digamma(Nk_out+tau);
//   NumericVector score(len_alpha+len_beta+1+1+1); // +1 for tau
//   long double zetak, L_k, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma, dLNk_dtau, dLk_dtau;
//   NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//   NumericMatrix cheese(len_alpha+len_beta+1+1+1,len_alpha+len_beta+1+1+1); // +1 for tau
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     zetak = zeta[k];
//
//     A1k = 0;
//     B1k = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2, eta, zeta
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//     dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//     ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//     for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//     ddeltaki_detaki[i] = 1;
//       for(j = (0); j < (len_beta); j++){
//         ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//       }
//       ddeltaki_sigma[i] = 0;
//     }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dalpha = 0*dLk_dalpha;
//     dLk_dgamma = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     dLk_dtau = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//     // likelihood and derivative contributions form size
//     //L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//     L_Nk = exp(lgam_Nktau[k]-logNk_factorial[k]-lgam_tau+tau*log(tau)+Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-(tau+Nk_out[k])*log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       for(j = (0); j < (len_alpha); j++){
//         //dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//         dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       }
//       //dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       //dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       //
//       dLNk_dtau = L_Nk*(DG_Nktau[k]-DG_tau+log(tau)+1-log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau)-(tau+Nk_out[k])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau) );
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//       L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k += GH_w[q]*L_Nk*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_alpha); j++){
//         dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
//       }
//       dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
//       for(j = (0); j < (len_beta); j++){
//         dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
//       }
//       dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma); //multiply by sigma because max wrt logsigma
//       dLk_dtau += GH_w[q]*tau*L_Yk*dLNk_dtau; //multiply by tau because max wrt logtau
//
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//     score[j] = dLk_dalpha[j]/(L_k+1e-323);
//     }
//     score[len_alpha] = dLk_dgamma/(L_k+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+1] = dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_alpha+1+len_beta] = dLk_dsigma/(L_k+1e-323);
//     score[len_alpha+1+len_beta+1] = dLk_dtau/(L_k+1e-323);
//
//
//     cheese += square_vect(score);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return cheese;
// }
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // tau as well!
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_JMMICS_NegBin_slopes(NumericVector gamma,
//                                            NumericVector sigma,
//                                            double tau,
//                                            double DG_tau,
//                                            NumericVector Omega,
//                                            NumericVector Delta,
//                                            NumericVector zeta,
//                                            NumericVector eta,
//                                            NumericVector rand_indic,
//                                            NumericVector Nk,
//                                            NumericVector Nk_actual,
//                                            NumericMatrix Zk,
//                                            NumericVector Yki,
//                                            NumericMatrix Xki,
//                                            NumericVector IDk,
//                                            NumericVector Wki,
//                                            int minNk,
//                                            bool condSize,
//                                            bool condOut,
//                                            NumericVector GH_w,
//                                            NumericVector GH_z,
//                                            int len_alpha,
//                                            int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector lgam_Nktau = lgamma(Nk_out+tau);
//   double lgam_tau = lgamma(tau);
//   NumericVector DG_Nktau = digamma(Nk_out+tau);
//   NumericVector score(len_alpha+len_beta+2+2+1); // +1 for tau //2 gammas and 2 sigmas
//   long double zetak, L_k, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_dgamma, dLNk_dsigma, dLYk_dsigma, dLk_dgamma, dLk_dsigma, dLNk_dtau, dLk_dtau;
//   NumericVector dLNk_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//   NumericMatrix cheese(len_alpha+len_beta+2+2+1,len_alpha+len_beta+2+2+1); // +1 for tau //2 gammas and 2 sigmas
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     zetak = zeta[k];
//
//     A1k = 0;
//     B1k = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2, eta, zeta
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//     dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dalpha = 0*dLk_dalpha;
//     dLk_dgamma = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     dLk_dtau = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions form size
//       //L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//       L_Nk = exp(lgam_Nktau[k]-logNk_factorial[k]-lgam_tau+tau*log(tau)+Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-(tau+Nk_out[k])*log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       for(j = (0); j < (len_alpha); j++){
//         //dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//         dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       }
//       //dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       //dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-(tau+Nk_out[k])*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//       //
//       dLNk_dtau = L_Nk*(DG_Nktau[k]-DG_tau+log(tau)+1-log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau)-(tau+Nk_out[k])/(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau) );
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k += GH_w[q]*L_Nk*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_alpha); j++){
//         dLk_dalpha[j] += GH_w[q]*L_Yk*dLNk_dalpha[j];
//       }
//       dLk_dgamma += GH_w[q]*L_Yk*dLNk_dgamma;
//       for(j = (0); j < (len_beta); j++){
//         dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk;
//       }
//       dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk+L_Yk*dLNk_dsigma); //multiply by sigma because max wrt logsigma
//       dLk_dtau += GH_w[q]*tau*L_Yk*dLNk_dtau; //multiply by tau because max wrt logtau
//
//
//
//     }
//
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//       score[j] = dLk_dalpha[j]/(L_k+1e-323);
//     }
//     score[len_alpha] = (1-rand_indic[k])*dLk_dgamma/(L_k+1e-323);
//     score[len_alpha+1] = rand_indic[k]*dLk_dgamma/(L_k+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+2] = dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_alpha+2+len_beta] = (1-rand_indic[k])*dLk_dsigma/(L_k+1e-323);
//     score[len_alpha+2+len_beta+1] = rand_indic[k]*dLk_dsigma/(L_k+1e-323);
//     score[len_alpha+2+len_beta+2] = dLk_dtau/(L_k+1e-323);
//
//
//
//     cheese += square_vect(score);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return cheese;
// }
//
//
// // [[Rcpp::export]]
// NumericVector compute_EB_JMMICS_NegBin(NumericVector gamma,
//                                        NumericVector sigma,
//                                       double tau,
//                                       NumericVector Omega,
//                                       NumericVector Delta,
//                                       NumericVector Nk,
//                                       NumericVector Nk_actual,
//                                       NumericVector Yki,
//                                       NumericVector IDk,
//                                       NumericVector Wki,
//                                       int minNk,
//                                       NumericVector GH_w,
//                                       NumericVector GH_z) {
//
//   int maxK = Omega.size(); // number of clusters
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i;
//   double L_Nk, L_Yk, L_k, zL_k;
//   NumericVector pred(maxK);
//   NumericVector Nk_out = Nk - minNk; // adjust for Poisson + minimum
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     L_k = 0; // kth denominator
//     zL_k = 0; // kth numerator
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       //L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//       L_Nk = exp(lgamma(Nk_out[k]+tau)-logNk_factorial[k]-lgamma(tau)+tau*log(tau)+Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-(tau+Nk_out[k])*log(exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])+tau));
//
//       L_Yk = 1; // contribution of Yk to likelihood
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//
//       }
//
//       L_k += GH_w[q]*L_Nk*L_Yk; // denominator
//       zL_k += GH_w[q]*L_Nk*L_Yk*sigma[k]*GH_z[q]; // numerator
//     }
//
//     cumNk += Nk_actual[k];
//
//     pred[k] = zL_k/(L_k+1e-323);
//   }
//
//   return pred;
// }
//
//
//
//
//
//
//
//
//
//
//
//
// //
// // ZIP MODEL FUNCTIONS //
// //
//
// // [[Rcpp::export]]
// double compute_logLik_JMMICS_ZIP(NumericVector gamma,
//                                  NumericVector sigma,
//                                  NumericVector Xi,
//                                  NumericVector Omega,
//                                  NumericVector Delta,
//                                  NumericVector Nk,
//                                  NumericVector Nk_actual,
//                                  NumericVector Nk0,
//                                  NumericVector Yki,
//                                  NumericVector IDk,
//                                  NumericVector Wki,
//                                  int minNk,
//                                  NumericVector GH_w,
//                                  NumericVector GH_z) {
//
//     int maxK = max(IDk); // number of clusters
//     int nquad = GH_w.size(); // number of quadrature points
//     int cumNk = 0; // to only subset outcomes from cluster k
//     long double logL_ZIP = 0; // log Likelihood (objective function to be maximized)
//     int k, q, i;
//     long double L_k_ZIP, L_Nk, L_Nk_ZIP, L_Yk, lam;
//     NumericVector Nk_out = Nk - minNk; // adjust for Poisson + minimum
//     NumericVector logNk_factorial = lgamma(Nk_out+1);
//
//     for(k = 0; k < maxK; k++) { // loop over clusters
//
//     L_k_ZIP = 0; // kth contribution to the (non-log) likelihood
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//     // its really logLNk
//     L_Nk =  Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q]) -exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k];   // contribution of Nk to likelihood
//     lam =  exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//     //L_Nk =  ::Rf_dpois(Nk_out[k],lam,1);
//     L_Nk_ZIP = log((Nk0[k])*expit(Xi[k])+(1-expit(Xi[k]))*exp(L_Nk)); // unfortunate to need to exponentiate and relog due to the sum  // can make it an if statement if necessary
//
//     L_Yk = 1; // contribution of Yk to likelihood
//
//     for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//     L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//
//     }
//
//     L_k_ZIP +=  exp(log(GH_w[q])+(L_Nk_ZIP)+log(L_Yk));
//     }
//
//     cumNk += Nk_actual[k];
//
//     logL_ZIP += log(L_k_ZIP+1e-323);
//     }
//
//     return logL_ZIP;
// }
//
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericVector compute_score_JMMICS_ZIP(NumericVector gamma,
//                                        NumericVector sigma,
//                                        NumericVector Xi,
//                                        NumericVector Omega,
//                                        NumericVector Delta,
//                                        NumericVector theta,
//                                        NumericVector zeta,
//                                        NumericVector eta,
//                                        NumericVector Nk,
//                                        NumericVector Nk_actual,
//                                        NumericVector Nk0,
//                                        NumericMatrix Z0k,
//                                        NumericMatrix Zk,
//                                        NumericVector Yki,
//                                        NumericMatrix Xki,
//                                        NumericVector IDk,
//                                        NumericVector Wki,
//                                        int minNk,
//                                        bool condSize,
//                                        bool condOut,
//                                        NumericVector GH_w,
//                                        NumericVector GH_z,
//                                        int len_epsilon,
//                                        int len_alpha,
//                                        int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector score(len_alpha+1+len_beta+1+len_epsilon);
//   long double zetak, thetak, L_k_ZIP, L_Nk_ZIP, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector dXik_depsilon(len_epsilon);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_ZIP_dgamma, dLNk_dgamma, dLNk_ZIP_dsigma, dLNk_dsigma, dLYk_dsigma, dLk_ZIP_dgamma, dLk_ZIP_dsigma;
//   NumericVector dLNk_ZIP_dalpha(len_alpha), dLNk_dalpha(len_alpha), dLk_ZIP_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLNk_ZIP_depsilon(len_epsilon), dLk_ZIP_depsilon(len_epsilon);
//   NumericVector dLYk_dbeta(len_beta), dLk_ZIP_dbeta(len_beta);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     zetak = zeta[k];
//     thetak = theta[k];
//
//     A1k = 0;
//     B1k = 0;
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2 A0 B0, eta, zeta, theta
//
//     // derivative chunks for zero model
//
//     for(j = (0); j < (len_epsilon); j++){
//       dXik_depsilon[j] =  Z0k(k,j);
//     }
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//       dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k_ZIP = 0;
//     dLk_ZIP_dalpha = 0*dLk_ZIP_dalpha;
//     dLk_ZIP_dgamma = 0;
//     dLk_ZIP_depsilon = 0*dLk_ZIP_depsilon;
//     dLk_ZIP_dbeta = 0*dLk_ZIP_dbeta;
//     dLk_ZIP_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from zero-inflated size model
//       L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of poisson to likelihood
//       L_Nk_ZIP = (Nk0[k])*expit(Xi[k])+(1-expit(Xi[k]))*(L_Nk); // full size model likelihood contribution
//       for(j = (0); j < (len_alpha); j++){
//         dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//         dLNk_ZIP_dalpha[j] =  (1-expit(Xi[k]))*dLNk_dalpha[j];
//       }
//       dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dgamma = (1-expit(Xi[k]))*dLNk_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLNk_ZIP_depsilon[j] =  dXik_depsilon[j]*expit(Xi[k])*(1-expit(Xi[k]))*(Nk0[k]-L_Nk);
//       }
//
//       dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dsigma = (1-expit(Xi[k]))*dLNk_dsigma;
//
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k_ZIP += GH_w[q]*L_Nk_ZIP*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_alpha); j++){
//         dLk_ZIP_dalpha[j] += GH_w[q]*L_Yk*dLNk_ZIP_dalpha[j];
//       }
//       dLk_ZIP_dgamma += GH_w[q]*L_Yk*dLNk_ZIP_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLk_ZIP_depsilon[j] += GH_w[q]*L_Yk*dLNk_ZIP_depsilon[j];
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLk_ZIP_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk_ZIP;
//       }
//       dLk_ZIP_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk_ZIP+L_Yk*dLNk_ZIP_dsigma);
//
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//       score[j] += dLk_ZIP_dalpha[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha] += dLk_ZIP_dgamma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+1] += dLk_ZIP_dbeta[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha+1+len_beta] += dLk_ZIP_dsigma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_epsilon); j++){
//       score[j+len_alpha+1+len_beta+1] += dLk_ZIP_depsilon[j]/(L_k_ZIP+1e-323);
//     }
//
//
//     cumNk += Nk_actual[k];
//   }
//
//   return score;
// }
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericVector compute_score_JMMICS_ZIP_slopes(NumericVector gamma,
//                                               NumericVector sigma,
//                                               NumericVector Xi,
//                                               NumericVector Omega,
//                                               NumericVector Delta,
//                                               NumericVector theta,
//                                               NumericVector zeta,
//                                               NumericVector eta,
//                                               NumericVector rand_indic,
//                                               NumericVector Nk,
//                                               NumericVector Nk_actual,
//                                               NumericVector Nk0,
//                                               NumericMatrix Z0k,
//                                               NumericMatrix Zk,
//                                               NumericVector Yki,
//                                               NumericMatrix Xki,
//                                               NumericVector IDk,
//                                               NumericVector Wki,
//                                               int minNk,
//                                               bool condSize,
//                                               bool condOut,
//                                               NumericVector GH_w,
//                                               NumericVector GH_z,
//                                               int len_epsilon,
//                                               int len_alpha,
//                                               int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector score(len_alpha+2+len_beta+2+len_epsilon); //2 gammas and 2 sigmas
//   long double zetak, thetak, L_k_ZIP, L_Nk_ZIP, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector dXik_depsilon(len_epsilon);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_ZIP_dgamma, dLNk_dgamma, dLNk_ZIP_dsigma, dLNk_dsigma, dLYk_dsigma, dLk_ZIP_dgamma, dLk_ZIP_dsigma;
//   NumericVector dLNk_ZIP_dalpha(len_alpha), dLNk_dalpha(len_alpha), dLk_ZIP_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLNk_ZIP_depsilon(len_epsilon), dLk_ZIP_depsilon(len_epsilon);
//   NumericVector dLYk_dbeta(len_beta), dLk_ZIP_dbeta(len_beta);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     zetak = zeta[k];
//     thetak = theta[k];
//
//     A1k = 0;
//     B1k = 0;
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2 A0 B0, eta, zeta, theta
//
//     // derivative chunks for zero model
//
//     for(j = (0); j < (len_epsilon); j++){
//       dXik_depsilon[j] =  Z0k(k,j);
//     }
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//       dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k_ZIP = 0;
//     dLk_ZIP_dalpha = 0*dLk_ZIP_dalpha;
//     dLk_ZIP_dgamma = 0;
//     dLk_ZIP_depsilon = 0*dLk_ZIP_depsilon;
//     dLk_ZIP_dbeta = 0*dLk_ZIP_dbeta;
//     dLk_ZIP_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from zero-inflated size model
//       L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of poisson to likelihood
//       L_Nk_ZIP = (Nk0[k])*expit(Xi[k])+(1-expit(Xi[k]))*(L_Nk); // full size model likelihood contribution
//       for(j = (0); j < (len_alpha); j++){
//         dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//         dLNk_ZIP_dalpha[j] =  (1-expit(Xi[k]))*dLNk_dalpha[j];
//       }
//       dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dgamma = (1-expit(Xi[k]))*dLNk_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLNk_ZIP_depsilon[j] =  dXik_depsilon[j]*expit(Xi[k])*(1-expit(Xi[k]))*(Nk0[k]-L_Nk);
//       }
//
//       dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dsigma = (1-expit(Xi[k]))*dLNk_dsigma;
//
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k_ZIP += GH_w[q]*L_Nk_ZIP*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_alpha); j++){
//         dLk_ZIP_dalpha[j] += GH_w[q]*L_Yk*dLNk_ZIP_dalpha[j];
//       }
//       dLk_ZIP_dgamma += GH_w[q]*L_Yk*dLNk_ZIP_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLk_ZIP_depsilon[j] += GH_w[q]*L_Yk*dLNk_ZIP_depsilon[j];
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLk_ZIP_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk_ZIP;
//       }
//       dLk_ZIP_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk_ZIP+L_Yk*dLNk_ZIP_dsigma);
//
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//       score[j] += dLk_ZIP_dalpha[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha] += (1-rand_indic[k])*dLk_ZIP_dgamma/(L_k_ZIP+1e-323);
//     score[len_alpha+1] += rand_indic[k]*dLk_ZIP_dgamma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+2] += dLk_ZIP_dbeta[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha+2+len_beta] += (1-rand_indic[k])*dLk_ZIP_dsigma/(L_k_ZIP+1e-323);
//     score[len_alpha+2+len_beta+1] += rand_indic[k]*dLk_ZIP_dsigma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_epsilon); j++){
//       score[j+len_alpha+2+len_beta+2] += dLk_ZIP_depsilon[j]/(L_k_ZIP+1e-323);
//     }
//
//
//
//     cumNk += Nk_actual[k];
//   }
//
//   return score;
// }
//
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_JMMICS_ZIP(NumericVector gamma,
//                                         NumericVector sigma,
//                                         NumericVector Xi,
//                                         NumericVector Omega,
//                                         NumericVector Delta,
//                                         NumericVector theta,
//                                         NumericVector zeta,
//                                         NumericVector eta,
//                                         NumericVector Nk,
//                                         NumericVector Nk_actual,
//                                         NumericVector Nk0,
//                                         NumericMatrix Z0k,
//                                         NumericMatrix Zk,
//                                         NumericVector Yki,
//                                         NumericMatrix Xki,
//                                         NumericVector IDk,
//                                         NumericVector Wki,
//                                         int minNk,
//                                         bool condSize,
//                                         bool condOut,
//                                         NumericVector GH_w,
//                                         NumericVector GH_z,
//                                         int len_epsilon,
//                                         int len_alpha,
//                                         int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector score(len_alpha+1+len_beta+1+len_epsilon);
//   long double zetak, thetak, L_k_ZIP, L_Nk_ZIP, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector dXik_depsilon(len_epsilon);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_ZIP_dgamma, dLNk_dgamma, dLNk_ZIP_dsigma, dLNk_dsigma, dLYk_dsigma, dLk_ZIP_dgamma, dLk_ZIP_dsigma;
//   NumericVector dLNk_ZIP_dalpha(len_alpha), dLNk_dalpha(len_alpha), dLk_ZIP_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLNk_ZIP_depsilon(len_epsilon), dLk_ZIP_depsilon(len_epsilon);
//   NumericVector dLYk_dbeta(len_beta), dLk_ZIP_dbeta(len_beta);
//   NumericMatrix cheese(len_alpha+1+len_beta+1+len_epsilon,len_alpha+1+len_beta+1+len_epsilon);
//
//
//
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//
//     zetak = zeta[k];
//     thetak = theta[k];
//
//     A1k = 0;
//     B1k = 0;
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2 A0 B0, eta, zeta, theta
//
//     // derivative chunks for zero model
//
//     for(j = (0); j < (len_epsilon); j++){
//       dXik_depsilon[j] =  Z0k(k,j);
//     }
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//       dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k_ZIP = 0;
//     dLk_ZIP_dalpha = 0*dLk_ZIP_dalpha;
//     dLk_ZIP_dgamma = 0;
//     dLk_ZIP_depsilon = 0*dLk_ZIP_depsilon;
//     dLk_ZIP_dbeta = 0*dLk_ZIP_dbeta;
//     dLk_ZIP_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from zero-inflated size model
//       L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of poisson to likelihood
//       L_Nk_ZIP = (Nk0[k])*expit(Xi[k])+(1-expit(Xi[k]))*(L_Nk); // full size model likelihood contribution
//       for(j = (0); j < (len_alpha); j++){
//         dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//         dLNk_ZIP_dalpha[j] =  (1-expit(Xi[k]))*dLNk_dalpha[j];
//       }
//       dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dgamma = (1-expit(Xi[k]))*dLNk_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLNk_ZIP_depsilon[j] =  dXik_depsilon[j]*expit(Xi[k])*(1-expit(Xi[k]))*(Nk0[k]-L_Nk);
//       }
//
//       dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dsigma = (1-expit(Xi[k]))*dLNk_dsigma;
//
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k_ZIP += GH_w[q]*L_Nk_ZIP*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_alpha); j++){
//         dLk_ZIP_dalpha[j] += GH_w[q]*L_Yk*dLNk_ZIP_dalpha[j];
//       }
//       dLk_ZIP_dgamma += GH_w[q]*L_Yk*dLNk_ZIP_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLk_ZIP_depsilon[j] += GH_w[q]*L_Yk*dLNk_ZIP_depsilon[j];
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLk_ZIP_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk_ZIP;
//       }
//       dLk_ZIP_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk_ZIP+L_Yk*dLNk_ZIP_dsigma);
//
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//       score[j] = dLk_ZIP_dalpha[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha] = dLk_ZIP_dgamma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+1] = dLk_ZIP_dbeta[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha+1+len_beta] = dLk_ZIP_dsigma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_epsilon); j++){
//       score[j+len_alpha+1+len_beta+1] = dLk_ZIP_depsilon[j]/(L_k_ZIP+1e-323);
//     }
//
//
//
//     cheese += square_vect(score);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return cheese;
// }
//
//
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_JMMICS_ZIP_slopes(NumericVector gamma,
//                                                NumericVector sigma,
//                                                NumericVector Xi,
//                                                NumericVector Omega,
//                                                NumericVector Delta,
//                                                NumericVector theta,
//                                                NumericVector zeta,
//                                                NumericVector eta,
//                                                NumericVector rand_indic,
//                                                NumericVector Nk,
//                                                NumericVector Nk_actual,
//                                                NumericVector Nk0,
//                                                NumericMatrix Z0k,
//                                                NumericMatrix Zk,
//                                                NumericVector Yki,
//                                                NumericMatrix Xki,
//                                                NumericVector IDk,
//                                                NumericVector Wki,
//                                                int minNk,
//                                                bool condSize,
//                                                bool condOut,
//                                                NumericVector GH_w,
//                                                NumericVector GH_z,
//                                                int len_epsilon,
//                                                int len_alpha,
//                                                int len_beta) {
//
//   int maxK = Omega.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector Nk_out = Nk - minNk;
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//   NumericVector score(len_alpha+2+len_beta+2+len_epsilon); //2 gammas and 2 sigmas
//   long double zetak, thetak, L_k_ZIP, L_Nk_ZIP, L_Nk, L_Yk, A1k, B1k;
//   NumericVector A2ki(N), B2ki(N);
//   long double dOmegak_dzetak, dOmegak_dgamma, dOmegak_dsigma;
//   NumericVector dOmegak_dalpha(len_alpha);
//   NumericVector dXik_depsilon(len_epsilon);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLNk_ZIP_dgamma, dLNk_dgamma, dLNk_ZIP_dsigma, dLNk_dsigma, dLYk_dsigma, dLk_ZIP_dgamma, dLk_ZIP_dsigma;
//   NumericVector dLNk_ZIP_dalpha(len_alpha), dLNk_dalpha(len_alpha), dLk_ZIP_dalpha(len_alpha), dLk_dalpha(len_alpha);
//   NumericVector dLNk_ZIP_depsilon(len_epsilon), dLk_ZIP_depsilon(len_epsilon);
//   NumericVector dLYk_dbeta(len_beta), dLk_ZIP_dbeta(len_beta);
//   NumericMatrix cheese(len_alpha+2+len_beta+2+len_epsilon,len_alpha+2+len_beta+2+len_epsilon); //2 gammas and 2 sigmas
//
//
//
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//
//     zetak = zeta[k];
//     thetak = theta[k];
//
//     A1k = 0;
//     B1k = 0;
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       A1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]);
//
//       B1k += GH_w[q]*exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])*GH_z[q];
//
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A1 B1 A2 B2 A0 B0, eta, zeta, theta
//
//     // derivative chunks for zero model
//
//     for(j = (0); j < (len_epsilon); j++){
//       dXik_depsilon[j] =  Z0k(k,j);
//     }
//
//     // derivative chunks for size model
//     if (condSize == FALSE) {
//       dOmegak_dgamma = -sigma[k]*B1k/(A1k+1e-323);
//       dOmegak_dsigma = -gamma[k]*B1k/(A1k+1e-323);
//       dOmegak_dzetak = exp(zetak)/(A1k+1e-323);
//     } else{  // much simpler if conditional model is used
//       dOmegak_dgamma = 0;
//       dOmegak_dsigma = 0;
//       dOmegak_dzetak = 1;
//     }
//     for(j = (0); j < (len_alpha); j++){
//       dOmegak_dalpha[j] =  dOmegak_dzetak*Zk(k,j);
//     }
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k_ZIP = 0;
//     dLk_ZIP_dalpha = 0*dLk_ZIP_dalpha;
//     dLk_ZIP_dgamma = 0;
//     dLk_ZIP_depsilon = 0*dLk_ZIP_depsilon;
//     dLk_ZIP_dbeta = 0*dLk_ZIP_dbeta;
//     dLk_ZIP_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from zero-inflated size model
//       L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of poisson to likelihood
//       L_Nk_ZIP = (Nk0[k])*expit(Xi[k])+(1-expit(Xi[k]))*(L_Nk); // full size model likelihood contribution
//       for(j = (0); j < (len_alpha); j++){
//         dLNk_dalpha[j] =  L_Nk*dOmegak_dalpha[j]*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//         dLNk_ZIP_dalpha[j] =  (1-expit(Xi[k]))*dLNk_dalpha[j];
//       }
//       dLNk_dgamma = L_Nk*(dOmegak_dgamma+sigma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dgamma = (1-expit(Xi[k]))*dLNk_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLNk_ZIP_depsilon[j] =  dXik_depsilon[j]*expit(Xi[k])*(1-expit(Xi[k]))*(Nk0[k]-L_Nk);
//       }
//
//       dLNk_dsigma = L_Nk*(dOmegak_dsigma+gamma[k]*GH_z[q])*(Nk_out[k]-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q]));
//       dLNk_ZIP_dsigma = (1-expit(Xi[k]))*dLNk_dsigma;
//
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k_ZIP += GH_w[q]*L_Nk_ZIP*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_alpha); j++){
//         dLk_ZIP_dalpha[j] += GH_w[q]*L_Yk*dLNk_ZIP_dalpha[j];
//       }
//       dLk_ZIP_dgamma += GH_w[q]*L_Yk*dLNk_ZIP_dgamma;
//       for(j = (0); j < (len_epsilon); j++){
//         dLk_ZIP_depsilon[j] += GH_w[q]*L_Yk*dLNk_ZIP_depsilon[j];
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLk_ZIP_dbeta[j] += GH_w[q]*dLYk_dbeta[j]*L_Nk_ZIP;
//       }
//       dLk_ZIP_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma*L_Nk_ZIP+L_Yk*dLNk_ZIP_dsigma);
//
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_alpha); j++){ // summing over the length of the score (ie num of params)
//       score[j] = dLk_ZIP_dalpha[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha] = (1-rand_indic[k])*dLk_ZIP_dgamma/(L_k_ZIP+1e-323);
//     score[len_alpha+1] = rand_indic[k]*dLk_ZIP_dgamma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_beta); j++){
//       score[j+len_alpha+2] = dLk_ZIP_dbeta[j]/(L_k_ZIP+1e-323);
//     }
//     score[len_alpha+2+len_beta] = (1-rand_indic[k])*dLk_ZIP_dsigma/(L_k_ZIP+1e-323);
//     score[len_alpha+2+len_beta+1] = rand_indic[k]*dLk_ZIP_dsigma/(L_k_ZIP+1e-323);
//     for(j = 0; j < (len_epsilon); j++){
//       score[j+len_alpha+2+len_beta+2] = dLk_ZIP_depsilon[j]/(L_k_ZIP+1e-323);
//     }
//
//
//
//     cheese += square_vect(score);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return cheese;
// }
//
//
//
//
//
//
// // [[Rcpp::export]]
// NumericVector compute_EB_JMMICS_ZIP(NumericVector gamma,
//                                     NumericVector sigma,
//                                     NumericVector Xi,
//                                     NumericVector Omega,
//                                     NumericVector Delta,
//                                     NumericVector Nk,
//                                     NumericVector Nk_actual,
//                                     NumericVector Nk0,
//                                     NumericVector Yki,
//                                     NumericVector IDk,
//                                     NumericVector Wki,
//                                     int minNk,
//                                     NumericVector GH_w,
//                                     NumericVector GH_z) {
//
//   int maxK = Omega.size(); // number of clusters
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i;
//   double L_Nk, L_Nk_ZIP, L_Yk, L_k_ZIP, zL_k_ZIP;
//   NumericVector pred(maxK);
//   NumericVector Nk_out = Nk - minNk; // adjust for Poisson + minimum
//   NumericVector logNk_factorial = lgamma(Nk_out+1);
//
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     L_k_ZIP = 0; // kth denominator
//     zL_k_ZIP = 0; // kth numerator
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       L_Nk = exp(Nk_out[k]*(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-exp(Omega[k]+gamma[k]*sigma[k]*GH_z[q])-logNk_factorial[k]);   // contribution of Nk to likelihood
//       L_Nk_ZIP = ((Nk0[k])*expit(Xi[k])+(1-expit(Xi[k]))*(L_Nk)); // unfortunate to need to exponentiate and relog due to the sum  // can make it an if statement if necessary
//
//       L_Yk = 1; // contribution of Yk to likelihood
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//
//       }
//
//       L_k_ZIP += GH_w[q]*L_Nk_ZIP*L_Yk; // denominator
//       zL_k_ZIP += GH_w[q]*L_Nk_ZIP*L_Yk*sigma[k]*GH_z[q]; // numerator
//     }
//
//     cumNk += Nk_actual[k];
//
//     pred[k] = zL_k_ZIP/(L_k_ZIP+1e-323);
//   }
//
//   return pred;
// }
//
//
//
//
//
//
//
//
// //
// // OUTCOME-ONLY MODEL //
// //
//
// // [[Rcpp::export]]
// double compute_logLik_GLMM(NumericVector sigma,
//                            NumericVector Delta,
//                            NumericVector Nk_actual,
//                            NumericVector Yki,
//                            NumericVector IDk,
//                            NumericVector Wki,
//                            NumericVector GH_w,
//                            NumericVector GH_z) {
//
//   int maxK = Nk_actual.size(); // number of clusters
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   long double logL = 0; // log Likelihood (objective function to be maximized)
//   int k, q, i;
//   long double L_k, L_Yk;
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     L_k = 0; // kth contribution to the (non-log) likelihood
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       L_Yk = 1; // contribution of Yk to likelihood
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//
//       }
//
//       L_k +=  exp(log(GH_w[q])+log(L_Yk));//exp(log(GH_w[q])+L_Nk+log(L_Yk));
//     }
//
//     cumNk += Nk_actual[k];
//
//     logL += log(L_k+1e-323);
//   }
//
//   return logL;
// }
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericVector compute_score_GLMM(NumericVector sigma,
//                                  NumericVector Delta,
//                                  NumericVector eta,
//                                  NumericVector Nk_actual,
//                                  NumericVector Yki,
//                                  NumericMatrix Xki,
//                                  NumericVector IDk,
//                                  NumericVector Wki,
//                                  bool condOut,
//                                  NumericVector GH_w,
//                                  NumericVector GH_z,
//                                  int len_beta) {
//
//   int maxK = Nk_actual.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector score(len_beta+1);
//   long double L_k, L_Yk;
//   NumericVector A2ki(N), B2ki(N);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLYk_dsigma, dLk_dsigma;
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A2 B2, eta
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k += GH_w[q]*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_beta); j++){
//         dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j];
//       }
//       dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma);
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_beta); j++){
//       score[j] += dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_beta] += dLk_dsigma/(L_k+1e-323);
//
//
//     cumNk += Nk_actual[k];
//   }
//
//   return score;
// }
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericVector compute_score_GLMM_slopes(NumericVector sigma,
//                                         NumericVector Delta,
//                                         NumericVector eta,
//                                         NumericVector rand_indic,
//                                         NumericVector Nk_actual,
//                                         NumericVector Yki,
//                                         NumericMatrix Xki,
//                                         NumericVector IDk,
//                                         NumericVector Wki,
//                                         bool condOut,
//                                         NumericVector GH_w,
//                                         NumericVector GH_z,
//                                         int len_beta) {
//
//   int maxK = Nk_actual.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector score(len_beta+2); //2 gammas and 2 sigmas
//   long double L_k, L_Yk;
//   NumericVector A2ki(N), B2ki(N);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLYk_dsigma, dLk_dsigma;
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A2 B2, eta
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k += GH_w[q]*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_beta); j++){
//         dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j];
//       }
//       dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma);
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_beta); j++){
//       score[j] += dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_beta] += (1-rand_indic[k])*dLk_dsigma/(L_k+1e-323);
//     score[len_beta+1] += rand_indic[k]*dLk_dsigma/(L_k+1e-323);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return score;
// }
//
//
//
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_GLMM(NumericVector sigma,
//                                   NumericVector Delta,
//                                   NumericVector eta,
//                                   NumericVector Nk_actual,
//                                   NumericVector Yki,
//                                   NumericMatrix Xki,
//                                   NumericVector IDk,
//                                   NumericVector Wki,
//                                   bool condOut,
//                                   NumericVector GH_w,
//                                   NumericVector GH_z,
//                                   int len_beta) {
//
//   int maxK = Nk_actual.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector score(len_beta+1);
//   long double L_k, L_Yk;
//   NumericVector A2ki(N), B2ki(N);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLYk_dsigma, dLk_dsigma;
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//   NumericMatrix cheese(len_beta+1,len_beta+1);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//        for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A2 B2, eta
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k += GH_w[q]*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_beta); j++){
//         dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j];
//       }
//       dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma);
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_beta); j++){
//       score[j] = dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_beta] = dLk_dsigma/(L_k+1e-323);
//
//
//     cheese += square_vect(score);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return cheese;
// }
//
// // note that the maximization occurs with respect to logsigma, so we multiply the Us by sigma for chain rule
// // [[Rcpp::export]]
// NumericMatrix compute_cheese_GLMM_slopes(NumericVector sigma,
//                                          NumericVector Delta,
//                                          NumericVector eta,
//                                          NumericVector rand_indic,
//                                          NumericVector Nk_actual,
//                                          NumericVector Yki,
//                                          NumericMatrix Xki,
//                                          NumericVector IDk,
//                                          NumericVector Wki,
//                                          bool condOut,
//                                          NumericVector GH_w,
//                                          NumericVector GH_z,
//                                          int len_beta) {
//
//   int maxK = Nk_actual.size(); // number of clusters
//   int N = Delta.size(); // number of obs
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i, j;
//   NumericVector score(len_beta+2); // 2 gammas and 2 sigmas
//   long double L_k, L_Yk;
//   NumericVector A2ki(N), B2ki(N);
//   NumericVector ddeltaki_detaki(N), ddeltaki_sigma(N);
//   NumericMatrix ddeltaki_dbeta(N,len_beta);
//   long double dLYk_dsigma, dLk_dsigma;
//   NumericVector dLYk_dbeta(len_beta), dLk_dbeta(len_beta);
//   NumericMatrix cheese(len_beta+2,len_beta+2); //2 sigmas
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         A2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]));
//
//         B2ki[i] += GH_w[q]*expit(Delta[i]+sigma[k]*GH_z[q])*(1-expit(Delta[i]+sigma[k]*GH_z[q]))*GH_z[q];
//
//       }
//
//     }
//
//     // at this point we have A2 B2, eta
//
//     // derivative chunks for outcome model
//     if (condOut == FALSE) {
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = expit(eta[i])*(1-expit(eta[i]))/(A2ki[i]+1e-323);
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = -B2ki[i]/(A2ki[i]+1e-323);
//       }
//     } else { // much simpler if conditional model is used
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//         ddeltaki_detaki[i] = 1;
//         for(j = (0); j < (len_beta); j++){
//           ddeltaki_dbeta(i,j) =  ddeltaki_detaki[i]*Xki(i,j);
//         }
//         ddeltaki_sigma[i] = 0;
//       }
//     }
//
//     // at this point we have deriv chunks that will be needed for the next round of integration
//     L_k = 0;
//     dLk_dbeta = 0*dLk_dbeta;
//     dLk_dsigma = 0;
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       // likelihood and derivative contributions from outcome
//       L_Yk = 1; // contribution of Yk to likelihood
//       dLYk_dsigma = 0;
//       dLYk_dbeta = 0*dLYk_dbeta;
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//         for(j = (0); j < (len_beta); j++){
//           dLYk_dbeta[j] += Wki[i]*ddeltaki_dbeta(i,j)*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//         }
//         dLYk_dsigma += Wki[i]*(ddeltaki_sigma[i]+GH_z[q])*(Yki[i]-expit(Delta[i]+sigma[k]*GH_z[q]));
//       }
//       for(j = (0); j < (len_beta); j++){
//         dLYk_dbeta[j] =  L_Yk*dLYk_dbeta[j];  // scaling it by L_Yk as required
//       }
//       dLYk_dsigma = L_Yk*dLYk_dsigma; // scaling it by L_Yk as required
//
//
//
//       L_k += GH_w[q]*L_Yk;
//       // compute main derivs here
//       for(j = (0); j < (len_beta); j++){
//         dLk_dbeta[j] += GH_w[q]*dLYk_dbeta[j];
//       }
//       dLk_dsigma += GH_w[q]*sigma[k]*(dLYk_dsigma);
//
//
//     }
//
//     //then scale all by L_k and sum over k
//     for(j = 0; j < (len_beta); j++){
//       score[j] = dLk_dbeta[j]/(L_k+1e-323);
//     }
//     score[len_beta] = (1-rand_indic[k])*dLk_dsigma/(L_k+1e-323);
//     score[len_beta+1] = rand_indic[k]*dLk_dsigma/(L_k+1e-323);
//
//
//
//     cheese += square_vect(score);
//
//     cumNk += Nk_actual[k];
//   }
//
//   return cheese;
// }
//
// // [[Rcpp::export]]
// NumericVector compute_EB_GLMM(NumericVector sigma,
//                               NumericVector Delta,
//                               NumericVector Nk_actual,
//                               NumericVector Yki,
//                               NumericVector IDk,
//                               NumericVector Wki,
//                               NumericVector GH_w,
//                               NumericVector GH_z) {
//
//   int maxK = Nk_actual.size(); // number of clusters
//   int nquad = GH_w.size(); // number of quadrature points
//   int cumNk = 0; // to only subset outcomes from cluster k
//   int k, q, i;
//   double L_Yk, L_k, zL_k;
//   NumericVector pred(maxK);
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     L_k = 0; // kth denominator
//     zL_k = 0; // kth numerator
//
//     for(q = 0; q < nquad; q++){ // loop over GH zeros
//
//       L_Yk = 1; // contribution of Yk to likelihood
//
//       for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//         L_Yk *=  pow(exp(Yki[i]*(Delta[i]+sigma[k]*GH_z[q]))/(1+exp(Delta[i]+sigma[k]*GH_z[q])), Wki[i]);
//
//       }
//
//       L_k += GH_w[q]*L_Yk; // denominator
//       zL_k += GH_w[q]*L_Yk*sigma[k]*GH_z[q]; // numerator
//     }
//
//     cumNk += Nk_actual[k];
//
//     pred[k] = zL_k/(L_k+1e-323);
//   }
//
//   return pred;
// }
//










// // compute naive EEE for C Code
// // [[Rcpp::export]]
// NumericVector compute_naiveEEE(NumericVector muki,
//                                NumericVector mu0ki,
//                                NumericVector mu1ki,
//                                NumericVector muXki,
//                                NumericMatrix Zki,
//                                NumericMatrix Zxki,
//                                int len_beta,
//                                int len_alphax,
//                                NumericVector Xki,
//                                NumericVector Nk_actual,
//                                NumericVector Yki,
//                                NumericVector IDk,
//                                NumericVector WEIGHTki,
//                                NumericVector VALIDki) {
//
//   int maxK = max(IDk); // number of clusters
//   int cumNk = 0; // to only subset outcomes from cluster k
//   long double EEEYki = 0; // EEE temporary variable
//   long double EEEXki = 0; // EEE temporary variable
//   NumericVector EEE(len_beta+len_alphax); // EEE vector
//   int k, i, p;
//
//
//   for(k = 0; k < maxK; k++) { // loop over clusters
//
//     for(i = (0+cumNk); i < (Nk_actual[k]+cumNk); i++){ // loop over outcomes
//
//       EEEYki = 0; // reset temporary variable
//       EEEXki = 0;
//
//       if(VALIDki[i]>0){ // validation sample
//         // EEEk += WEIGHTki[i]*(Yki[i]-expit(muki[i]))*(Xki[i]*expit(muXki[i]) + (1-Xki[i])*(1-expit(muXki[i]))); // dont think we need the P(X|) chunk of this since we can just use the comlete data function
//         EEEYki += WEIGHTki[i]*(Yki[i]-expit(muki[i])); //complete data EEE // *(Xki[i]*expit(muXki[i]) + (1-Xki[i])*(1-expit(muXki[i])));
//         EEEXki += (Xki[i]-expit(muXki[i])); // NO WEIGHTS? (N is in the model)
//
//         // multiply by covariate and add to vector
//         EEE[0] += Zki(i,0)*EEEYki; //intercept
//         EEE[1] += Xki[i]*EEEYki; // x term
//         for(p = 2; p < (len_beta); p++){ // loop over parameters
//           EEE[p] += Zki(i,p-1)*EEEYki;  //-1 not -2 because we've done the intercept, but X isnt in Zki
//         }
//         for(p = len_beta; p < (len_beta+len_alphax); p++){ // loop over parameters
//           EEE[p] += Zxki(i,p-len_beta)*EEEXki;  // - len beta since were only looking at Zxki now
//         }
//
//
//       } else { // marginalize in main sample
//         EEEYki += WEIGHTki[i]*(Yki[i]-expit(mu0ki[i]))*(1-expit(muXki[i])); // X=0
//         EEEYki += WEIGHTki[i]*(Yki[i]-expit(mu1ki[i]))*(expit(muXki[i]));  // X=1
//
//         // these cancel to 0 for main sample
//         // EEEXki += (0-expit(muXki[i]))*(1-expit(muXki[i])); // X=0 // NO WEIGHTS? (N is in the model)
//         // EEEXki += (1-expit(muXki[i]))*(expit(muXki[i]));  // X=1
//
//         // multiply by covariate and add to vector
//         EEE[0] += Zki(i,0)*EEEYki; //intercept
//         EEE[1] += WEIGHTki[i]*(Yki[i]-expit(mu1ki[i]))*(expit(muXki[i])); // second term is 0 //Xki[i]*EEEYki; // x term
//         for(p = 2; p < (len_beta); p++){ // loop over parameters
//           EEE[p] += Zki(i,p-1)*EEEYki;  //-1 not -2 because we've done the intercept, but X isnt in Zki
//         }
//         // for(p = len_beta; p < (len_beta+len_alphax); p++){ // loop over parameters
//         //   // These will all be 0---could comment out if need be
//         //   EEE[p] += Zxki(i,p)*EEEXki;  //-1 not -2 because we've done the intercept, but X isnt in Zki
//         // }
//       }
//
//
//
//
//     }
//
//   }
//
//
//   return EEE;
//
// }








