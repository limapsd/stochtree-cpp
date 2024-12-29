/*! Copyright (c) 2024 stochtree authors. All rights reserved. */
#ifndef STOCHTREE_ALPHA_SAMPLER_H_
#define STOCHTREE_ALPHA_SAMPLER_H_

#include <random>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace StochTree{
class AlphaSampler{
 public:
    AlphaSampler(){}
   ~AlphaSampler(){}
  //  double Sample(std::vector<double>& lpv, double& alpha, double a, double b, double rho, std::mt19937 & gen){
  //   size_t p = lpv.size();
  //   double sumlpv = 0.;
    
  //   std::vector<double> lambda_g (1000,0.);
  //   std::vector<double> theta_g (1000,0.);
  //   std::vector<double> loglikes (1000,0.);
  //   for(size_t j=0;j<p;j++) sumlpv+=lpv[j];
    
  //   for(size_t k=0;k<1000;k++){
  //     lambda_g[k]=(double)(k+1)/1001.;
    
  //     theta_g[k]=(lambda_g[k]*rho)/(1.-lambda_g[k]);
    
  //     double theta_log_lik=lgamma(theta_g[k])-(double)p*lgamma(theta_g[k]/(double)p)+(theta_g[k]/(double)p)*sumlpv;
    
  //     double beta_log_prior=(a-1.)*log(lambda_g[k])+(b-1.)*log(1.-lambda_g[k]);
    
  //     loglikes[k]=theta_log_lik+beta_log_prior;      
  //   }
  //   double log_sum_exp = logsumexp(loglikes);
    
  //   for(size_t k =0; k<1000; k++){
  //       loglikes[k] = exp(loglikes[k] - log_sum_exp);
  //   }
  //   int index = sample_class(loglikes,gen);
  //   double r_alpha = theta_g[index];

  //   return r_alpha;
  //  }
  std::tuple< std::vector<double>, double, double> Sample(std::vector<double>& lpv, double a, double b, double rho, std::mt19937 & gen){
    size_t p = lpv.size();
    double sumlpv = 0.;
    
    std::vector<double> lambda_g (1000,0.);
    std::vector<double> theta_g (1000,0.);
    std::vector<double> loglikes (1000,0.);
    for(size_t j=0;j<p;j++) sumlpv+=lpv[j];
    
    for(size_t k=0;k<1000;k++){
      lambda_g[k]=(double)(k+1)/1001.;
    
      theta_g[k]=(lambda_g[k]*rho)/(1.-lambda_g[k]);
    
      double theta_log_lik=lgamma(theta_g[k])-(double)p*lgamma(theta_g[k]/(double)p)+(theta_g[k]/(double)p)*sumlpv;
    
      double beta_log_prior=(a-1.)*log(lambda_g[k])+(b-1.)*log(1.-lambda_g[k]);
    
      loglikes[k]=theta_log_lik+beta_log_prior;      
    }
    double log_sum_exp = logsumexp(loglikes);
    
    for(size_t k =0; k<1000; k++){
        loglikes[k] = exp(loglikes[k] - log_sum_exp);
    }
    int index = sample_class(loglikes,gen);
    double r_alpha = theta_g[index];

    return std::make_tuple(loglikes, r_alpha, log_sum_exp);
   }
   
 private:
    double logsumexp(std::vector<double>& v){
            double sm = 0.;
            double mx = *std::max_element(v.begin(), v.end());
            for(int i = 0 ;i< v.size(); i++){
                sm += exp(v[i] - mx);
            }
            return mx + log(sm);
        }
    int sample_class(std::vector<double>& probs, std::mt19937& gen) {
      
      std::uniform_real_distribution<double>unif_rand(0,1);
      double U = unif_rand(gen);
      double foo = 0.0;
      int K = probs.size();

    // Sample
    for(int k = 0; k < K; k++) {
        foo += probs[k];
        if(U < foo) {
        return(k);
        }
    }
    return K - 1;
    }
      
};

}

#endif //STOCHTREE_ALPHA_SAMPLER_H_