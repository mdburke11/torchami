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

 TamiBase::g_prod_t construct_example2();
 TamiBase::ami_vars construct_ext_example2(TamiBase& tami, int ebatch_size, int fbatch_size);

 TamiBase::ami_vars construct_4ord_ext_multipole_example(TamiBase& tami, int ebatch_size, int fbatch_size);
 
 TamiBase::g_prod_t construct_multipole_example();

 TamiBase::ami_vars construct_ext_example6(TamiBase& tami, int ebatch_size, int fbatch_size);
 TamiBase::g_prod_t construct_example6();


 TamiBase::ami_vars construct_ext_example1_bose(TamiBase& tami, int ebatch_size, int fbatch_size);
 TamiBase::g_prod_t construct_example1_bose();
 
 TamiBase::g_prod_t construct_example_Y();
 TamiBase::ami_vars construct_ext_example_Y(TamiBase& tami, int ebatch_size, int fbatch_size);

 TamiBase::g_prod_t construct_example_J();
 TamiBase::ami_vars construct_ext_example_J(TamiBase& tami, int ebatch_size, int fbatch_size);

 TamiBase::g_prod_t construct_example6();
 TamiBase::ami_vars construct_ext_example6(TamiBase& tami, int ebatch_size, int fbatch_size);

 void default_example();
 void example2();
 void example1_bose();
 void example4();
 void example6();
 void example9();
 void graph_library_example();
 void renorm_PT_graph_example();

// TODO: clean these up - maybe put together
 std::string format_r1_tensor(const at::Tensor&);
 std::string format_r2_tensor(const at::Tensor&);
