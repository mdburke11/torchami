#include "examples.hpp"

TamiBase::g_prod_t construct_example2() {

  TamiBase::g_prod_t g;

  // Setting up G array
  // defining alpha's /// std::vector<int>

  TamiBase::alpha_t alpha_1 = {1, 0, 0};
  TamiBase::alpha_t alpha_2 = {0, 1, 0};
  TamiBase::alpha_t alpha_3 = {-1, 1, 1};

  // defining epsilon's
  TamiBase::epsilon_t epsilon_1 = {1, 0, 0};
  TamiBase::epsilon_t epsilon_2 = {0, 1, 0};
  TamiBase::epsilon_t epsilon_3 = {0, 0, 1};

  TamiBase::g_struct g1(epsilon_1, alpha_1);
  TamiBase::g_struct g2(epsilon_2, alpha_2);
  TamiBase::g_struct g3(epsilon_3, alpha_3);

  // TamiBase::g_prod_t R0={g1,g2,g3};

  // OR
  TamiBase::g_prod_t R0;
  R0.push_back(g1);
  R0.push_back(g2);
  R0.push_back(g3);

  return R0;
}

TamiBase::ami_vars construct_ext_example2(TamiBase &tami, int ebatch_size,
                                          int fbatch_size) {

  std::vector<TamiBase::complex_double> energy_vec = {
      -4, 0.1, -1}; // template energy vector that will be copied batch_size
                    // times in the tensor
  int energy_size = energy_vec.size();
  TamiBase::energy_t energy =
      8.0 * at::rand({ebatch_size, energy_size}, tami.options) - 4.0;
  std::vector<at::Tensor> freq_vecs = {};
  for (int i = 0; i < fbatch_size; ++i) {
    freq_vecs.push_back(at::tensor(
        {TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, M_PI * (2 * i + 1))},
        tami.options));
  }

  TamiBase::frequency_t frequency = at::vstack(freq_vecs);

  double BETA = 1.0;
  TamiBase::ami_vars external(energy, frequency, BETA);

  return external;
}

TamiBase::ami_vars construct_ext_example_J(TamiBase &tami, int ebatch_size,
                                           int fbatch_size) {

  std::vector<TamiBase::complex_double> energy_vec = {
      -4.64, 1.02,  1.04,  1.05, 1.06, 1.07, 1.08, 1.09, 1.11,
      1.23,  -4.43, -4.52, 1.5,  1.6,  1.7,  1.8,  1.9};
  int energy_size = energy_vec.size();
  TamiBase::energy_t energy =
      8.0 * at::rand({ebatch_size, energy_size}, tami.options) - 4.0;
  std::vector<at::Tensor> freq_vecs = {};
  for (int i = 0; i < fbatch_size; ++i) {
    freq_vecs.push_back(at::tensor(
        {TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, M_PI * (2 * i + 1))},
        tami.options));
  }

  TamiBase::frequency_t frequency = at::vstack(freq_vecs);

  double BETA = 1.0;
  TamiBase::ami_vars external(energy, frequency, BETA);

  return external;
}

TamiBase::g_prod_t construct_example_J() {

  TamiBase::g_prod_t g;

  // Setting up G array
  // defining alpha's

  TamiBase::alpha_t alpha_1 = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_2 = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_3 = {1, 1, 0, 0, 0, 0, 0, 0, 0, -1};
  TamiBase::alpha_t alpha_4 = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_5 = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_6 = {1, 1, 1, -1, 0, 0, 0, 0, 0, -1};
  TamiBase::alpha_t alpha_7 = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_8 = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
  TamiBase::alpha_t alpha_9 = {0, -1, 0, 0, 1, 0, 0, 1, 0, 0};
  TamiBase::alpha_t alpha_10 = {0, 0, 1, -1, 1, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_11 = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_12 = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
  TamiBase::alpha_t alpha_13 = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
  TamiBase::alpha_t alpha_14 = {0, 0, -1, 1, -1, 1, 1, 0, 0, 0};
  TamiBase::alpha_t alpha_15 = {0, -1, -1, 1, 0, 0, 0, 1, 0, 1};
  TamiBase::alpha_t alpha_16 = {0, 0, 1, -1, 1, -1, 0, 0, 1, 0};
  TamiBase::alpha_t alpha_17 = {0, -1, 0, 0, 1, -1, 0, 1, 1, 0};

  // defining epsilon's
  TamiBase::epsilon_t epsilon_1 = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_2 = {0, 1, 0, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_3 = {0, 0, 1, 0, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_4 = {0, 0, 0, 1, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_5 = {0, 0, 0, 0, 1, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_6 = {0, 0, 0, 0, 0, 1, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_7 = {0, 0, 0, 0, 0, 0, 1, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_8 = {0, 0, 0, 0, 0, 0, 0, 1, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_9 = {0, 0, 0, 0, 0, 0, 0, 0, 1,
                                   0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_10 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    1, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_11 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 1, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_12 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 1, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_13 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 1, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_14 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 1, 0, 0, 0};
  TamiBase::epsilon_t epsilon_15 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 1, 0, 0};
  TamiBase::epsilon_t epsilon_16 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 1, 0};
  TamiBase::epsilon_t epsilon_17 = {0, 0, 0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0, 1};

  TamiBase::g_struct g1(epsilon_1, alpha_1);
  TamiBase::g_struct g2(epsilon_2, alpha_2);
  TamiBase::g_struct g3(epsilon_3, alpha_3);
  TamiBase::g_struct g4(epsilon_4, alpha_4);
  TamiBase::g_struct g5(epsilon_5, alpha_5);
  TamiBase::g_struct g6(epsilon_6, alpha_6);
  TamiBase::g_struct g7(epsilon_7, alpha_7);
  TamiBase::g_struct g8(epsilon_8, alpha_8);
  TamiBase::g_struct g9(epsilon_9, alpha_9);
  TamiBase::g_struct g10(epsilon_10, alpha_10);
  TamiBase::g_struct g11(epsilon_11, alpha_11);
  TamiBase::g_struct g12(epsilon_12, alpha_12);
  TamiBase::g_struct g13(epsilon_13, alpha_13);
  TamiBase::g_struct g14(epsilon_14, alpha_14);
  TamiBase::g_struct g15(epsilon_15, alpha_15);
  TamiBase::g_struct g16(epsilon_16, alpha_16);
  TamiBase::g_struct g17(epsilon_17, alpha_17);

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
  R0.push_back(g12);
  R0.push_back(g13);
  R0.push_back(g14);
  R0.push_back(g15);
  R0.push_back(g16);
  R0.push_back(g17);

  return R0;
}

/// @brief construction function for a 4th order diagram to make ami_vars object
/// with ebatch_size and fbatch_size energies and frequencies
/// @param tami, ebatch_size, fbatch_size
/// @return TamiBase::ami_vars
TamiBase::ami_vars construct_4ord_ext_multipole_example(TamiBase &tami,
                                                        int ebatch_size,
                                                        int fbatch_size) {

  std::vector<TamiBase::complex_double> energy_vec = {
      1, 1.1, 1.2, 1.31, 1.4, 0.01, 0.1}; //{1,1.1,1.2,1.3,1.4,0.01, 0.1};
  int energy_size = energy_vec.size();

  TamiBase::energy_t energy =
      8.0 * at::rand({ebatch_size, energy_size}, tami.options) - 4.0;
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

TamiBase::g_prod_t construct_example1_bose() {

  TamiBase::g_prod_t g;

  // Setting up G array
  // defining alpha's /// std::vector<int>

  TamiBase::alpha_t alpha_1 = {1, 0};
  TamiBase::alpha_t alpha_2 = {1, 1};

  // defining epsilon's
  TamiBase::epsilon_t epsilon_1 = {1, 0};
  TamiBase::epsilon_t epsilon_2 = {0, 1};

  TamiBase::g_struct g1(epsilon_1, alpha_1);
  TamiBase::g_struct g2(epsilon_2, alpha_2);

  TamiBase::g_prod_t R0;

  R0.push_back(g1);
  R0.push_back(g2);

  return R0;
}

TamiBase::ami_vars construct_ext_example1_bose(TamiBase &tami, int ebatch_size,
                                               int fbatch_size) {

  std::vector<TamiBase::complex_double> energy_vec = {
      -4, 0.1}; // Residual from libami code
  int energy_size = energy_vec.size();

  TamiBase::energy_t energy =
      8.0 * at::rand({ebatch_size, energy_size}, tami.options) -
      4.0; // at::zeros({1, energy_size},
           // tami.options);//8.0*at::rand({ebatch_size,energy_size},
           // tami.options)-4.0;

  // energy[0][0]=-4;
  // energy[0][1]=0.1;

  std::vector<at::Tensor> freq_vecs = {};
  for (int i = 0; i < fbatch_size; ++i) {
    freq_vecs.push_back(
        at::tensor({TamiBase::complex_double(0, 0),
                    TamiBase::complex_double(0, M_PI * (2 * i + 1))},
                   tami.options));
  }

  TamiBase::frequency_t frequency = at::vstack(freq_vecs);
  double BETA = 1.0;
  TamiBase::ami_vars external(energy, frequency, BETA);

  return external;
}

TamiBase::ami_vars construct_ext_example_Y(TamiBase &tami, int ebatch_size,
                                           int fbatch_size) {

  std::vector<TamiBase::complex_double> energy_vec = {4, -1, -1, -1, -2, 2, 1};
  int energy_size = energy_vec.size();

  TamiBase::energy_t energy =
      8.0 * at::rand({ebatch_size, energy_size}, tami.options) - 4.0;

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

TamiBase::g_prod_t construct_example_Y() {

  TamiBase::g_prod_t g;

  // Setting up G array
  // defining alpha's

  TamiBase::alpha_t alpha_1 = {1, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_2 = {0, 1, 0, 0, 0};
  TamiBase::alpha_t alpha_3 = {0, 0, 1, 0, 0};
  TamiBase::alpha_t alpha_4 = {0, 0, 0, 1, 0};
  TamiBase::alpha_t alpha_5 = {-1, 1, 0, 0, 1};
  TamiBase::alpha_t alpha_6 = {0, 1, 1, -1, 0};
  TamiBase::alpha_t alpha_7 = {1, 0, 1, 0, -1};

  // defining epsilon's
  TamiBase::epsilon_t epsilon_1 = {1, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_2 = {0, 1, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_3 = {0, 0, 1, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_4 = {0, 0, 0, 1, 0, 0, 0};
  TamiBase::epsilon_t epsilon_5 = {0, 0, 0, 0, 1, 0, 0};
  TamiBase::epsilon_t epsilon_6 = {0, 0, 0, 0, 0, 1, 0};
  TamiBase::epsilon_t epsilon_7 = {0, 0, 0, 0, 0, 0, 1};

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

TamiBase::ami_vars construct_ext_example6(TamiBase &tami, int ebatch_size,
                                          int fbatch_size) {

  std::vector<TamiBase::complex_double> energy_vec = {
      1, 1.1, 1.2, 1.3, 1.4, 0, 0.1, 0.2, 0.3, 0.4, 0.5};
  int energy_size = energy_vec.size();

  TamiBase::energy_t energy =
      8.0 * at::rand({ebatch_size, energy_size}, tami.options) - 4.0;

  std::vector<at::Tensor> freq_vecs = {};
  for (int i = 0; i < fbatch_size; ++i) {
    freq_vecs.push_back(at::tensor(
        {TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, 0), TamiBase::complex_double(0, 0),
         TamiBase::complex_double(0, M_PI * (2 * i + 1))},
        tami.options));
  }

  TamiBase::frequency_t frequency = at::vstack(freq_vecs);

  double BETA = 1.0;
  TamiBase::ami_vars external(energy, frequency, BETA);

  return external;
}

TamiBase::g_prod_t construct_example6() {

  TamiBase::g_prod_t g;

  // Setting up G array
  // defining alpha's

  TamiBase::alpha_t alpha_1 = {1, 0, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_2 = {0, 1, 0, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_3 = {0, 0, 1, 0, 0, 0, 0};
  TamiBase::alpha_t alpha_4 = {0, 0, 0, 1, 0, 0, 0};
  TamiBase::alpha_t alpha_5 = {0, 0, 0, 0, 1, 0, 0};
  TamiBase::alpha_t alpha_6 = {0, 0, 0, 0, 0, 1, 0};
  TamiBase::alpha_t alpha_7 = {1, -1, 0, 0, 0, 1, 0};
  TamiBase::alpha_t alpha_8 = {-1, 1, 0, 0, 0, 0, 1};
  TamiBase::alpha_t alpha_9 = {-1, 1, 0, 0, 0, 0, 1};
  TamiBase::alpha_t alpha_10 = {-1, 1, -1, 1, 0, 0, 1};
  TamiBase::alpha_t alpha_11 = {1, 0, 0, 0, -1, 1, 0};

  // defining epsilon's
  TamiBase::epsilon_t epsilon_1 = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_2 = {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_3 = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_4 = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_5 = {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_6 = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_7 = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0};
  TamiBase::epsilon_t epsilon_8 = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
  TamiBase::epsilon_t epsilon_9 = {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
  TamiBase::epsilon_t epsilon_10 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
  TamiBase::epsilon_t epsilon_11 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

  TamiBase::g_struct g1(epsilon_1, alpha_1);
  TamiBase::g_struct g2(epsilon_2, alpha_2);
  TamiBase::g_struct g3(epsilon_3, alpha_3);
  TamiBase::g_struct g4(epsilon_4, alpha_4);
  TamiBase::g_struct g5(epsilon_5, alpha_5);
  TamiBase::g_struct g6(epsilon_6, alpha_6);
  TamiBase::g_struct g7(epsilon_7, alpha_7);
  TamiBase::g_struct g8(epsilon_8, alpha_8);
  TamiBase::g_struct g9(epsilon_9, alpha_9);
  TamiBase::g_struct g10(epsilon_10, alpha_10);
  TamiBase::g_struct g11(epsilon_11, alpha_11);

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

std::string format_r1_tensor(const at::Tensor &tens) {

  std::ostringstream str;
  str << std::setprecision(10);
  for (int i = 0; i < tens.size(0); ++i) {
    auto a = tens[i].item();
    str << a << std::endl;
  }
  str << "Stored on device: " << tens.device()
      << std::endl; // print the device the tensor is stored on
  str << std::endl;
  return str.str();
}

std::string format_r2_tensor(const at::Tensor &tens) {
  std::ostringstream str;
  for (int i = 0; i < tens.size(0); ++i) {
    for (int j = 0; j < tens.size(1); ++j) {
      auto a = tens[i][j].item();
      str << a << " ";
    }
    str << std::endl;
  }
  str << "Stored on device: " << tens.device()
      << std::endl; // print the device the tensor is stored on
  return str.str();
}