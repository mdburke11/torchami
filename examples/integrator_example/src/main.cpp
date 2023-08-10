#include "integration.hpp"

int main(){
    
    matsubara_freq_test();

    return 0;
}

at::Tensor epsilon_tight_binding(at::Tensor k){
    return -2 * at::cos(k).sum(/*dim=*/1, /*keepdim=*/1);
}

at::Tensor my_func(at::Tensor x){
    return at::cos(x) + c10::complex<double>(0, 1) * at::sin(x);
}


void matsubara_freq_test(){

    // init device and tamibase  and tamigraph obj
    at::Device myDev = at::kCUDA;
    TamiBase ami(myDev);
    TamiBase::graph_type graph_type = TamiBase::Sigma;
    int seed = 0;
    TamiGraph g(graph_type, seed);

    // init graph
    int ord=2; int group=0; int n=0;
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
    double reW = 0.0;
    double imW = (2*mat_freq + 1) * M_PI / beta;

    ext_vars evars(beta, mu, k, reW, imW);
    TamiBase::ami_vars avars = prep_ext(ord, evars, myDev);

    // helpers
    TamiBase::ft_terms ftout;
    double E_REG = 0;
    int N_INT = 2;
    TamiBase::ami_parms parms(N_INT, E_REG);
    ami.construct(N_INT, R0, ftout);

    /*
    // test ct function
    std::vector<complete_R0_graph> ct_graphs = read_ct_diagrams(g, second_ord_sigma.graph, ord, 2);

    for (auto ct_g : ct_graphs){
        print_R0(ct_g.R0);
    }
    */

    /*

    flat_mc_integrator::integ_domain domain = {{0, 1}};
    std::function fn = at::sin;


    int batch_size = 100000000;
    at::TensorOptions integOptions = at::TensorOptions().dtype(at::kComplexDouble).device(at::kCUDA);
    flat_mc_integrator mc(integOptions, 1000000);
    integration_result result = mc.integrate(fn, 1, batch_size, domain);

    std::cout << at::real(result.ans) << " " << at::real(result.error) << std::endl;

    */

    AMI_integrand integrand(ami, R0, avars, ftout, parms, epsilon_tight_binding, 1, evars);

    at::TensorOptions integOptions = at::TensorOptions().dtype(at::kComplexDouble).device(at::kCUDA);
    int max_batch_size = 100000;
    flat_mc_integrator mc(integOptions, max_batch_size);
    flat_mc_integrator::integ_domain domain = {{0, 2*M_PI}, {0, 2*M_PI}, {0, 2*M_PI},{0, 2*M_PI}};

    std::cout << std::setprecision(10) << std::endl;
    for (int n=0; n<20; ++n){
        evars.imW = (2*n + 1) * M_PI / evars.beta;
        integrand.update_ext_vars(evars);

        evars.print_ext_vars();
        integration_result result = mc.integrate(integrand, 4, 1000000, domain);
        std::cout << mat_freq << " " << evars.imW << " " << second_ord_sigma.prefactor * result.ans / std::pow(2* M_PI, 4) << " " << result.error/ std::pow(2* M_PI, 4) << std::endl;
    }

    

    //int batch_size = 100000000;
    //at::TensorOptions integOptions = at::TensorOptions().dtype(at::kComplexDouble).device(at::kCUDA);
    //flat_mc_integrator mc(integOptions, 1000000);

















}