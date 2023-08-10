#include "integration.hpp"


complete_graph read_diagram(TamiGraph tamig, std::string folder, int ord, int group, int n){
// reads in the specific graph described by ord, group, and n

    TamiGraph::gg_matrix_t ggm;

    tamig.read_ggmp(folder, ggm, ord, ord+1);
    tamig.print_ggm(ggm);

    std::cout << "After removing pairs:" << std::endl;
    tamig.print_ggm(ggm);
    tamig.ggm_label(ggm, 0);

    TamiGraph::graph_t diagram = ggm[ord][group].graph_vec[n];
    double prefactor = tamig.get_prefactor(diagram, ord);

    return complete_graph(diagram, prefactor);
}


std::vector<complete_R0_graph> read_ct_diagrams(TamiGraph tamig, TamiGraph::graph_t diagram, int order, int max_ct){

    std::vector<complete_R0_graph> results;

    TamiBase::g_prod_t R0; // temp to read in the vectors
    double pref;

    std::vector<TamiGraph::graph_t> graphs;

    tamig.generate_sigma_ct(diagram, graphs, max_ct);

    for (auto x: graphs){
        tamig.graph_to_R0(x, R0);
        double pref = tamig.get_prefactor(x, order);
        results.push_back({R0, pref});
    }
    return results;
}
TamiBase::ami_vars prep_ext(int ord, ext_vars evars, at::Device dev){

    int default_batchsize = 10;
    at::TensorOptions options = at::TensorOptions().dtype(at::kDouble).device(dev);

    TamiBase::energy_t energy = at::zeros({default_batchsize, ord + 1}, options);
    TamiBase::frequency_t frequency;
    for (int i=0; i < ord; ++i) {frequency.push_back((0.0, 0.0));}
    frequency.push_back((evars.reW, evars.imW));
    TamiBase::ami_vars external{energy, frequency, evars.beta};

    return external;

}

void print_R0(TamiBase::g_prod_t R0){

    std::cout << "Alpha: " << std::endl;

    for (TamiBase::g_struct x : R0){
        for (int a : x.alpha_){
            std::cout << a << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Epsilon: " << std::endl;

    for (TamiBase::g_struct x : R0){
        for (int e : x.eps_){
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

}
