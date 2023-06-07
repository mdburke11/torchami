#include "tami_base.hpp"

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

 TamiBase::g_prod_t construct_example2();
 TamiBase::ami_vars construct_ext_example2(TamiBase& tami);

 TamiBase::ami_vars construct_4ord_ext_multipole_example(TamiBase& tami);
 TamiBase::g_prod_t construct_multipole_example();
 
 TamiBase::ami_vars construct_ext_example1_bose(TamiBase& tami);
 TamiBase::g_prod_t construct_example1_bose();
 
 TamiBase::g_prod_t construct_example_Y();
 TamiBase::ami_vars construct_ext_example_Y(TamiBase& tami);

 TamiBase::g_prod_t construct_example_J();
 TamiBase::ami_vars construct_ext_example_J(TamiBase& tami);

 TamiBase::g_prod_t construct_example6();
 TamiBase::ami_vars construct_ext_example6(TamiBase& tami);

 void default_example();
 void example_1();
 void example2();
 void example1_bose();
 void example4();
 void example6();
 void example9();
 //void example_();

// TODO: clean these up - maybe put together
 std::string format_r1_tensor(const at::Tensor&);
 std::string format_r2_tensor(const at::Tensor&);
