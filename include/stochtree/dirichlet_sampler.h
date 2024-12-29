/*! Copyright (c) 2024 stochtree authors. All rights reserved. */
#ifndef STOCHTREE_DIRICHLET_SAMPLER_H_
#define STOCHTREE_DIRICHLET_SAMPLER_H_

#include <random>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

namespace StochTree{

class DirichletSampler{
 public:
    DirichletSampler() {}
    ~DirichletSampler() {}
    // std::vector<double>Sample(std::vector<double> alpha, std::mt19937 & gen){
    //     size_t n = alpha.size();
    //     std::vector<double> draw(n);
    //     std::vector<double> prob_vector(n);

    //     for(int j =0; j < n; j++){
    //         draw[j] = log_gamma_sampler(alpha[j], gen);
    //     }
    //     double lse = logsumexp(draw);
        
    //     for(int i = 0; i < n; i++){
    //         draw[i] -= lse; 
    //         // prob_vector[i] = exp(draw[i]);
    //     }
    //     return draw;
    // }
    std::tuple<double, std::vector<double> >Sample(std::vector<double> alpha, std::mt19937 & gen){
        size_t n = alpha.size();
        std::vector<double> draw_g(n);
        std::vector<double> draw_d(n);
        std::vector<double> prob_vector(n);

        for(int j =0; j < n; j++){
            draw_g[j] = log_gamma_sampler(alpha[j], gen);
        }
        double lse = logsumexp(draw_g);
        
        for(int i = 0; i < n; i++){
            draw_d[i] = draw_g[i] -  lse; 
            // prob_vector[i] = exp(draw[i]);
        }
        return std::make_tuple(lse, draw_d);
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
    double log_gamma_sampler(double shape, std::mt19937& gen){
        std::gamma_distribution<double> gamma_dist(shape + 1.0, 1.);
        std::uniform_real_distribution<double> unif_rand(0, 1.);

        double y = log(gamma_dist(gen));
        double z = log(unif_rand(gen))/shape;

        return y+z;
    }
};

} // namespace StochTree

#endif // STOCHTREE_DIRICHLET_SAMPLER_H_