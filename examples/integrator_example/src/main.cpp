#include "integration.hpp"

int main() {
  matsubara_freq_test();
  return 0;
}

at::Tensor epsilon_tight_binding(at::Tensor k) {
  return -2 * at::cos(k).sum(/*dim=*/1, /*keepdim=*/1);
}

at::Tensor my_func(at::Tensor x) {
  return at::cos(x) + c10::complex<double>(0, 1) * at::sin(x);
}

void matsubara_freq_test() {

  // output 2nd order diagram data to this file
  // Use the python script to plot it
  std::string ofname = "data.dat";
  std::ofstream out;
  out.open(ofname);

  // init device and tamibase  and tamigraph obj
  at::Device myDev = at::kCUDA;
  TamiBase ami(myDev);
  TamiBase::graph_type graph_type = TamiBase::Sigma;
  int seed = 0;
  TamiGraph g(graph_type, seed);

  // init graph
  int ord = 2;
  int group = 0;
  int n = 0;
  std::string folder = "../../ggm_examples/";

  complete_graph second_ord_sigma = read_diagram(g, folder, ord, group, n);
  TamiBase::g_prod_t R0;
  g.graph_to_R0(second_ord_sigma.graph, R0);

  print_R0(R0);

  // external vars
  int mat_freq = 1;
  double beta = 8.33;
  TamiBase::complex_double mu{0.1, 0};
  std::vector<double> k{M_PI, 0};
  std::cout << k[0] << " " << k[1] << k.size() << std::endl;

  int N_freq = 70;

  std ::vector<at ::Tensor> freq_vecs = {};
  for (int i = 0; i < N_freq; ++i) {
    at::Tensor row =
        at::full({1, ord + 1}, TamiBase::complex_double(0, 0), myDev);
    row[0][-1] = TamiBase::complex_double(0, M_PI * (2 * i + 1) / beta);
    freq_vecs.push_back(row);
  }
  TamiBase::frequency_t frequencies = at::vstack(freq_vecs);
  ext_vars evars(beta, mu, k);
  TamiBase::ami_vars avars = prep_ext(ord, evars, myDev);
  avars.frequency_ = frequencies;

  // helpers
  TamiBase::ft_terms ftout;
  double E_REG = 0;
  int N_INT = 2;
  TamiBase::ami_parms parms(N_INT, E_REG);
  ami.construct(N_INT, R0, ftout);

  AMI_integrand integrand(ami, R0, avars, ftout, parms, epsilon_tight_binding,
                          /*eval_realpart*/ 0, evars);

  at::TensorOptions integOptions =
      at::TensorOptions().dtype(at::kComplexDouble).device(at::kCUDA);
  int max_batch_size = 100000;
  flat_mc_integrator mc(integOptions, N_freq, max_batch_size);
  flat_mc_integrator::integ_domain domain = {
      {0, 2 * M_PI}, {0, 2 * M_PI}, {0, 2 * M_PI}, {0, 2 * M_PI}};

  std::cout << std::setprecision(10) << std::endl;

  auto t1 = std::chrono::high_resolution_clock::now();
  integration_result result = mc.integrate(integrand, 4, 1000000, domain);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> difft = t2 - t1;
  std::chrono::nanoseconds dt =
      std::chrono::duration_cast<std::chrono::nanoseconds>(difft);

  for (int n = 0; n < N_freq; ++n) {
    std::cout << n << " " << at::real(frequencies[n][-1]).item<double>() << " "
              << at::imag(frequencies[n][-1]).item<double>() << " "
              << second_ord_sigma.prefactor *
                     at::real(result.ans[0][n]).item<double>()
              << " " << at::real(result.error[0][n]).item<double>()
              << std::endl;
    out << n << " " << at::real(frequencies[n][-1]).item<double>() << " "
        << at::imag(frequencies[n][-1]).item<double>() << " "
        << second_ord_sigma.prefactor *
               at::real(result.ans[0][n]).item<double>()
        << " " << at::real(result.error[0][n]).item<double>() << "\n";
  }
  std::cout << std::endl
            << "time elapsed: " << dt.count() / std::pow(10, 9) << " seconds"
            << std::endl;

  out.close();
}