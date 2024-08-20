

#include "../src/tami_base_src/tami_base.hpp"
#include "../src/tami_graph_src/tami_graph.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

#include "gtest/gtest.h"
// include other headers if necessary   #include ../src/ami.hpps


#include <torch/extension.h>
#include <torch/python.h>
#include <torch/torch.h>

TEST(fermi_bose, fermi_test){

//(int m, double sigma, double beta, at::Tensor(TamiBase::complex_double E) )
// m: mth derivative
// sigma: +1 for fermion, -1 for boson
// beta: inverse template
// E: complex energy 

TamiBase obj;

TamiBase::complex_double from_func, from_analytic;

std::vector<at::Tensor> tensors(1, at::tensor({TamiBase::complex_double(2.0,0)},obj.options));
TamiBase::energy_t energy = at::vstack(tensors);


std::complex<double> E(2.0,0);
double beta=1.0;

from_func=obj.fermi_bose(0, 1.0,  beta, energy)[0][0].item<TamiBase::complex_double>();

from_analytic=1.0/(std::exp(beta*E)+1.0);


ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());

}	

TEST(fermi_bose, fermi_test_deriv1){
	
TamiBase obj;

TamiBase::complex_double from_func(0,0);
TamiBase::complex_double from_analytic(0,0);
TamiBase::complex_double E(2.0,0);


std::vector<at::Tensor> tensors(1, at::tensor({TamiBase::complex_double(2.0,0)},obj.options));
TamiBase::energy_t energy = at::vstack(tensors);
double beta=1.0;

from_func=obj.fermi_bose(1, 1.0,  beta, energy)[0][0].item<TamiBase::complex_double>();

from_analytic=-( beta*std::exp(E*beta) / ( std::pow( (std::exp(E*beta) +1.0), 2.0)) );

ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());	
	
}

TEST(fermi_bose, fermi_test_deriv2){
	
TamiBase obj;

TamiBase::complex_double from_func(0,0);
TamiBase::complex_double from_analytic(0,0);
TamiBase::complex_double E(2.0,0);


std::vector<at::Tensor> tensors(1, at::tensor({TamiBase::complex_double(2.0,0)},obj.options));
TamiBase::energy_t energy = at::vstack(tensors);
double beta=1.0;

from_func=obj.fermi_bose(2, 1.0,  beta, energy)[0][0].item<TamiBase::complex_double>();

from_analytic=2.0*(std::pow(beta,2)*std::exp(2.0*beta*E))/ ( std::pow( (std::exp(E*beta) +1.0), 3.0)) - std::pow(beta,2)*std::exp(beta*E)/ ( std::pow( (std::exp(E*beta) +1.0), 2.0));

ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());	
	
}

TEST(fermi_bose, bose_test_deriv1){
	
TamiBase obj;

// std::complex<double> from_func, from_analytic;
TamiBase::complex_double from_func(0,0);
TamiBase::complex_double from_analytic(0,0);
TamiBase::complex_double E(2.0,0);


std::vector<at::Tensor> tensors(1, at::tensor({TamiBase::complex_double(2.0,0)},obj.options));
TamiBase::energy_t energy = at::vstack(tensors);
double beta=1.0;

from_func=obj.fermi_bose(1, -1.0,  beta, energy)[0][0].item<TamiBase::complex_double>();

from_analytic=beta*std::exp(beta*E)/std::pow(std::exp(beta*E)-1.0,2);

ASSERT_DOUBLE_EQ(from_func.real(),from_analytic.real());	
	
}
