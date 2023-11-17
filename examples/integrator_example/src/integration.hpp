#pragma once
#include "tami_base.hpp"
#include "tami_graph.hpp"

#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <unistd.h>

#include <chrono>
#include <thread>

// This code does not take advantage of the simultaneous frequency evaluations!
// However this example code is still useful to see how one could implement a
// Monte carlo integration scheme.

// external functions

struct ext_vars{

    double beta;
    TamiBase::complex_double mu;
    std::vector<double> k;
    double reW, imW;

    ext_vars(){
        beta=0; mu=0; k={0, 0}; reW=0; imW=0;
    }

    ext_vars(double bta, TamiBase::complex_double mu_, std::vector<double> k_, double reW_, double imW_){
        beta = bta; mu = mu_, k = k_; reW = reW_, imW = imW_;
    }

    void print_ext_vars(){
        std::cout << beta << " " << mu << " ";
        for (auto q : k){std::cout << q << " ";}
        std::cout << reW << " " << imW << std::endl;
    }
};

at::Tensor epsilon_tight_binding(at::Tensor k);

// graph utilities

struct complete_graph{
    // struct to hold the graph and associated graph prefactor together in this example

    TamiGraph::graph_t graph; // graph
    double prefactor; // graph prefactor

    complete_graph(TamiGraph::graph_t graph_, double pref){
        graph = graph_;
        prefactor = pref;
    }
};

struct complete_R0_graph{
    // struct to hold the graph and associated graph prefactor together in this example

    TamiBase::g_prod_t R0; // graph
    double prefactor; // graph prefactor

    complete_R0_graph(TamiBase::g_prod_t R0_, double pref){
        R0 = R0_;
        prefactor = pref;
    }
};

complete_graph read_diagram(TamiGraph tamig, std::string folder, int ord, int group, int n);
std::vector<complete_R0_graph> read_ct_diagrams(TamiGraph tamig, TamiGraph::graph_t diagram, int order, int max_ct);
TamiBase::ami_vars prep_ext(int ord, ext_vars evars, at::Device dev);
void print_R0(TamiBase::g_prod_t R0);

// integrator functions

class integration_result{

    public:
        double sum;
        double sum2;
        int N;
        
        double ans;
        double error;

        integration_result(at::Tensor& sum_, at::Tensor& sum2_, int N_){
            sum=sum_.item<double>(); sum2=sum2_.item<double>(); N=N_;

            ans = sum / N;
            double avg_ans2 = sum2 / N;
            error = std::sqrt((avg_ans2 - std::pow(ans, 2)) / N);
        }

};

class flat_mc_integrator{

    public:
        at::TensorOptions options = at::TensorOptions().dtype(at::kDouble).device(at::kCPU);
        int max_batch;

        flat_mc_integrator(at::TensorOptions options_, int max_bs){
            options = options_; max_batch = max_bs;
        }

        typedef std::vector<std::vector<double>> integ_domain;


        at::Tensor prepInput(int length, integ_domain& domain);

        integration_result integrate(std::function<at::Tensor(at::Tensor)> fn, int dim, int N, integ_domain domain);

};


// Momentum integrand functions

class AMI_integrand{
    public:

        TamiBase tami;
        TamiBase::g_prod_t R0;
        TamiBase::ami_vars avars;
        TamiBase::ft_terms ft;
        TamiBase::ami_parms parms; 
        std::function<at::Tensor(at::Tensor)> eps;
        bool evalReal;
        ext_vars evars;

        int dim;
        at::Device device = at::kCPU; // default to cpu, but cuda after initializer
        at::TensorOptions options;
        at::Tensor full_alpha; // kron product of alpha and identity(size=dim of prob eg: 2D)
        int order; // number of integrals to evaluate


        AMI_integrand(TamiBase tami_, TamiBase::g_prod_t R0_, TamiBase::ami_vars avars_, 
                    TamiBase::ft_terms ft_, TamiBase::ami_parms parms_, 
                    std::function<at::Tensor(at::Tensor)> eps_, bool evalReal_, ext_vars evars_){
            tami = tami_;
            R0 = R0_;
            avars = avars_;
            ft = ft_;
            parms = parms_;
            eps = eps_;
            evalReal = evalReal_;
            evars = evars_;

            // insert external variables in the AMI objects
            this->update_ext_vars(evars);

            // useful things to extract
            dim = evars.k.size();
            device = tami.getDevice();

            options = at::TensorOptions().device(device).dtype(at::kDouble);
            at::Tensor I = at::eye(dim, options);

            // extract alpha matrix place in torch tensor
            std::vector<at::Tensor> alpha_vec;
            for (int i=0; i < R0.size(); ++i){
                std::vector<int> row = R0[i].alpha_;
                alpha_vec.push_back(at::tensor(row, options));
            }
            at::Tensor alpha = at::vstack(alpha_vec);

            full_alpha = at::transpose(at::kron(alpha, I), 0, 1);
            order = parms.N_INT_;
        }
    
        void update_ext_vars(ext_vars new_evars){
            // assuming theres only one frequency
            evars = new_evars;
            //avars.frequency_[0][-1] = TamiBase::complex_double(new_evars.reW, new_evars.imW); 
            avars.BETA_ = new_evars.beta; 
        }

        ext_vars get_ext_vars(){
            return evars;
        }

        void update_integrand(at::Tensor k){
            // Takes a new set of internal momenta (eg: k=rand(0, 2pi)) with external and gets
            // correct lin. comb then evaluates the dispersion for each propagator in the diagram
            // and inserts the values into the avars object

            // stack the columns of all the inputs together
            std::vector<at::Tensor> cols{k};
            for (double q : evars.k){cols.push_back(at::full({at::size(k, /*dim*/ 0), 1}, q, options));}
            at::Tensor K = at::hstack(cols);
            at::Tensor combined = at::matmul(K, full_alpha);

            std::vector<at::Tensor> e_cols;
            for (at::Tensor x : at::split(combined, /*size of split*/ dim, /*dim to split on*/ 1)){
                e_cols.push_back(eps(x));
            }

            avars.energy_ = evars.mu - at::hstack(e_cols);
        }

        at::Tensor operator()(at::Tensor x){
            this->update_integrand(x);
            at::Tensor value = at::reshape(tami.evaluate(parms, ft, avars), {-1});
            if (evalReal){
                return at::real(value);
            }
            return at::imag(value);
        }
};


// main functions
void matsubara_freq_test();
