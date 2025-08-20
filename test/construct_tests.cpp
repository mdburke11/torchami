
#include "tami_base.hpp"
#include "tami_graph.hpp"
#include <iomanip>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <random>

#include "gtest/gtest.h"

// #include <torch/extension.h>
// #include <torch/python.h>
#include <torch/torch.h>


TamiBase::g_prod_t construct_example_6();

TEST(construct_tests, term_deriv){
	
TamiBase ami;

//take_derivative_gprod(g_prod_t &g_prod, pole_struct pole, double start_sign, Ri_t &r_out, pole_array_t &poles, sign_t &signs);

TamiBase::alpha_t alpha_1={1,0,0,0,0,0,0};
TamiBase::alpha_t alpha_2={0,1,0,0,0,0,0};
TamiBase::alpha_t alpha_3={0,0,1,0,0,0,0};
	
TamiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0,0,0,0,0};

TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);
TamiBase::g_struct g3(epsilon_3,alpha_3);

// TamiBase::pole_struct p1(epsilon_1,alpha_1);



TamiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);

TamiBase::pole_array_t poles=ami.find_poles(0,R0);
TamiBase::pole_struct pole=poles[0];

TamiBase::term term;
term.g_list=R0;

TamiBase::terms outterms;

ami.take_term_derivative(term, pole, outterms );
	
ASSERT_EQ(outterms.size(),1);
ASSERT_EQ(outterms[0].g_list.size(),4);
	
}


TEST(construct_tests, term_deriv2){
	
TamiBase ami;

//take_derivative_gprod(g_prod_t &g_prod, pole_struct pole, double start_sign, Ri_t &r_out, pole_array_t &poles, sign_t &signs);

TamiBase::alpha_t alpha_1={1,0,0,0,0,0,0};
TamiBase::alpha_t alpha_2={1,1,0,0,0,0,0};
TamiBase::alpha_t alpha_3={0,0,1,0,0,0,0};
	
TamiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0,0,0,0,0};

TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);
TamiBase::g_struct g3(epsilon_3,alpha_3);




TamiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);


TamiBase::pole_array_t poles=ami.find_poles(0,R0);
TamiBase::pole_struct pole=poles[0];

TamiBase::term term;
term.g_list=R0;

TamiBase::terms outterms;

ami.take_term_derivative(term, pole, outterms );
ASSERT_EQ(outterms.size(),2);

ASSERT_EQ(outterms[0].g_list.size(),4);

}



TEST(construct_tests, update_g){
	
TamiBase ami;

TamiBase::alpha_t alpha_1={-1,1,0,0,0,0,1};	
TamiBase::epsilon_t eps_1={0,1,0,0,1,0,-1,0,0,0,0};

TamiBase::g_struct g1(eps_1,alpha_1);


TamiBase::g_prod_t R0=construct_example_6();
TamiBase::pole_array_t pole_list=ami.find_poles(0,R0);

TamiBase::g_struct updated=ami.update_G_pole(g1, pole_list[0]);

// print_g_struct_info(g1);

// print_g_struct_info(updated);

TamiBase::epsilon_t eps_result={1, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0};
TamiBase::alpha_t alpha_result={0, 1, 0, 0, 0, 0, 1};

TamiBase::pole_struct p1(eps_result,alpha_result);
TamiBase::pole_struct p2(updated.eps_,updated.alpha_);

ASSERT_EQ(ami.pole_equiv(p1, p2),true);
	
	
}

TEST(construct_tests, find_poles){
	

TamiBase ami;

// Problem setup (see ami_example.cpp)
TamiBase::g_prod_t R0=construct_example_6();

TamiBase::pole_array_t pole_list_0=ami.find_poles(0,R0);
TamiBase::pole_array_t pole_list_1=ami.find_poles(1,R0);
TamiBase::pole_array_t pole_list_2=ami.find_poles(2,R0);
TamiBase::pole_array_t pole_list_3=ami.find_poles(3,R0);
TamiBase::pole_array_t pole_list_6=ami.find_poles(6,R0);

std::cout<<"Expect warnings: { ";
TamiBase::pole_array_t pole_list_7=ami.find_poles(7,R0);
TamiBase::pole_array_t pole_list_12=ami.find_poles(12,R0);
std::cout<<"} Warnings caught"<<std::endl;
	
ASSERT_EQ(pole_list_0.size(),6);
ASSERT_EQ(pole_list_1.size(),6);
ASSERT_EQ(pole_list_2.size(),5);
ASSERT_EQ(pole_list_3.size(),4);
ASSERT_EQ(pole_list_6.size(),5);
ASSERT_EQ(pole_list_7.size(),0);
ASSERT_EQ(pole_list_12.size(),0);	
	
}






TamiBase::g_prod_t construct_example_6(){

TamiBase::g_prod_t g;


// Setting up G array
// defining alpha's

TamiBase::alpha_t alpha_1={1,0,0,0,0,0,0};
TamiBase::alpha_t alpha_2={0,1,0,0,0,0,0};
TamiBase::alpha_t alpha_3={0,0,1,0,0,0,0};
TamiBase::alpha_t alpha_4={0,0,0,1,0,0,0};
TamiBase::alpha_t alpha_5={0,0,0,0,1,0,0};
TamiBase::alpha_t alpha_6={0,0,0,0,0,1,0};
TamiBase::alpha_t alpha_7={-1,1,0,0,0,0,1};
TamiBase::alpha_t alpha_8={-1,1,1,0,0,0,1};
TamiBase::alpha_t alpha_9={-1,1,1,1,0,0,1};
TamiBase::alpha_t alpha_10={-1,1,1,1,1,0,1};
TamiBase::alpha_t alpha_11={-1,1,1,1,1,1,1};

//defining epsilon's
TamiBase::epsilon_t epsilon_1={1,0,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_2={0,1,0,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_3={0,0,1,0,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_4={0,0,0,1,0,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_5={0,0,0,0,1,0,0,0,0,0,0};
TamiBase::epsilon_t epsilon_6={0,0,0,0,0,1,0,0,0,0,0};
TamiBase::epsilon_t epsilon_7={0,0,0,0,0,0,1,0,0,0,0};
TamiBase::epsilon_t epsilon_8={0,0,0,0,0,0,0,1,0,0,0};
TamiBase::epsilon_t epsilon_9={0,0,0,0,0,0,0,0,1,0,0};
TamiBase::epsilon_t epsilon_10={0,0,0,0,0,0,0,0,0,1,0};
TamiBase::epsilon_t epsilon_11={0,0,0,0,0,0,0,0,0,0,1};


TamiBase::g_struct g1(epsilon_1,alpha_1);
TamiBase::g_struct g2(epsilon_2,alpha_2);
TamiBase::g_struct g3(epsilon_3,alpha_3);
TamiBase::g_struct g4(epsilon_4,alpha_4);
TamiBase::g_struct g5(epsilon_5,alpha_5);
TamiBase::g_struct g6(epsilon_6,alpha_6);
TamiBase::g_struct g7(epsilon_7,alpha_7);
TamiBase::g_struct g8(epsilon_8,alpha_8);
TamiBase::g_struct g9(epsilon_9,alpha_9);
TamiBase::g_struct g10(epsilon_10,alpha_10);
TamiBase::g_struct g11(epsilon_11,alpha_11);


TamiBase::g_prod_t R0;

R0.push_back(g1);
R0.push_back(g2);
R0.push_back(g3);
R0.push_back(g4);
R0.push_back(g5);
R0.push_back(g6);
R0.push_back(g7);
R0.push_back(g8);
R0.push_back(g9);
R0.push_back(g10);
R0.push_back(g11);





return R0;

}


