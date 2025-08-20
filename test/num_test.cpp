#include <gtest/gtest.h>

#include "tami_base.hpp"
#include "tami_graph.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

// #include <torch/extension.h>
// #include <torch/python.h>
#include <torch/torch.h>

at::Device get_device();

void example4(TamiBase::complex_double &value);
TamiBase::ami_vars construct_4ord_ext_multipole_example(TamiBase &tami,
                                                        int ebatch_size,
                                                        int fbatch_size);
TamiBase::g_prod_t construct_multipole_example();



// This is a basic test to check that gtest is functioning properly - due to the GLIBC flag this would fail and has no dependence on torchami
TEST(NumTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}


at::Device get_device(){

at::Device myDev(at::kCPU);
if (torch::cuda::is_available()) {at::Device myCudaDev(at::kCUDA); return myCudaDev;}

return myDev;

}


TEST(NumTest, o4_test){

at::Device myDev(get_device());

TamiBase::complex_double from_benchmark;
example4(from_benchmark);	

std::complex<double> from_func(-0.0009888959048396373,0.00075325682931293159);

std::cout<<from_benchmark<<" =? "<<from_func<<std::endl;
	
ASSERT_NEAR(from_func.real(),from_benchmark.real(),1e-11);
ASSERT_NEAR(from_func.imag(),from_benchmark.imag(),1e-11);
}



// 

// void example4(std::complex<double> &value);
// at::Device get_device();


// 





void example4(TamiBase::complex_double &value){

//START example	
// class instance 
int fbatch_size = 1;
int ebatch_size = 100;

at::Device myDev(get_device());

std::cout<<"Test using device "<< myDev<<std::endl;

TamiBase ami(myDev);

// std::cout<<"Device set"<<std::endl;

// Problem setup
TamiBase::g_prod_t R0=construct_multipole_example();
TamiBase::ami_vars avars=construct_4ord_ext_multipole_example(ami, ebatch_size, fbatch_size);

	
// Integration/Evaluation parameters
double E_REG = 0; // Numerical regulator for small energies.  If inf/nan
				// results try E_REG=1e-8
int N_INT = 4;    // Number of Matsubara sums to perform
TamiBase::ami_parms test_amiparms(N_INT, E_REG);


// simplified storage type
TamiBase::ft_terms amiterms;

// Construct solution for problem defined in R0
ami.construct(N_INT, R0, amiterms);

// Evaluate solution
at::Tensor term_val = ami.evaluate(test_amiparms, amiterms,
									avars);

// std::cout<< term_val[0][0]<<std::endl;
value=term_val[0][0].item<TamiBase::complex_double>();


return ;
	
}



TamiBase::ami_vars construct_4ord_ext_multipole_example(TamiBase &tami,
                                                        int ebatch_size,
                                                        int fbatch_size) {

//   at::tensor energy_vec = at::tensor({1, 1.1, 1.2, 1.31, 1.4, 0.01, 0.1}, tami.options); 
//   int energy_size = energy_vec.size();

// 	std::cout<< energy_vec<<std::endl;


  std::vector<at::Tensor> tensors(ebatch_size, at::tensor({TamiBase::complex_double(1,0), TamiBase::complex_double(1.1,0), TamiBase::complex_double(1.2,0), TamiBase::complex_double(1.31,0), TamiBase::complex_double(1.4,0), TamiBase::complex_double(0.01,0), TamiBase::complex_double(0.1,0)},tami.options));

  TamiBase::energy_t energy = at::vstack(tensors);



  std::vector<at::Tensor> freq_vecs = {};
  for (int i = 0; i < fbatch_size; ++i) {
    freq_vecs.push_back(at::tensor(
        {TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, M_PI * (2 * i + 1))},
        tami.options));
  }

  TamiBase::frequency_t frequency = at::vstack(freq_vecs);

  double BETA = 1.0;
  TamiBase::ami_vars external(energy, frequency, BETA);

  return external;
}

TamiBase::g_prod_t construct_multipole_example() {

  TamiBase::g_prod_t g;

  // Setting up G array
  // defining alpha's
  // TamiBase::alpha_t alpha_1={1,0,0,1,-1};
  // TamiBase::alpha_t alpha_2={0,0,0,1,0};
  // TamiBase::alpha_t alpha_3={1,0,0,0,0};
  // TamiBase::alpha_t alpha_4={1,0,0,0,0};
  // TamiBase::alpha_t alpha_5={0,1,0,0,0};
  // TamiBase::alpha_t alpha_6={-1,1,1,0,0};
  // TamiBase::alpha_t alpha_7={0,0,1,0,0};

  TamiBase::alpha_t alpha_1 = {0, 0, 1, 1, -1};
  TamiBase::alpha_t alpha_2 = {0, 0, 0, 1, 0};
  TamiBase::alpha_t alpha_3 = {0, 0, 1, 0, 0};
  TamiBase::alpha_t alpha_4 = {0, 0, 1, 0, 0};
  TamiBase::alpha_t alpha_5 = {0, 1, 0, 0, 0};
  TamiBase::alpha_t alpha_6 = {1, 1, -1, 0, 0};
  TamiBase::alpha_t alpha_7 = {1, 0, 0, 0, 0};

  // defining epsilon's
  TamiBase::epsilon_t epsilon_1 = {0, 0, 0, 0, 0, 1, 0};
  TamiBase::epsilon_t epsilon_2 = {0, 0, 0, 0, 1, 0, 0};
  TamiBase::epsilon_t epsilon_3 = {1, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_4 = {1, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_5 = {0, 1, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_6 = {0, 0, 0, 1, 0, 0, 0};
  TamiBase::epsilon_t epsilon_7 = {0, 0, 1, 0, 0, 0, 0};

  TamiBase::g_struct g1(epsilon_1, alpha_1);
  TamiBase::g_struct g2(epsilon_2, alpha_2);
  TamiBase::g_struct g3(epsilon_3, alpha_3);
  TamiBase::g_struct g4(epsilon_4, alpha_4);
  TamiBase::g_struct g5(epsilon_5, alpha_5);
  TamiBase::g_struct g6(epsilon_6, alpha_6);
  TamiBase::g_struct g7(epsilon_7, alpha_7);

  TamiBase::g_prod_t R0;

  R0.push_back(g1);
  R0.push_back(g2);
  R0.push_back(g3);
  R0.push_back(g4);
  R0.push_back(g5);
  R0.push_back(g6);
  R0.push_back(g7);

  return R0;
}

