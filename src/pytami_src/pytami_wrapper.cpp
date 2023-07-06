#include "../src/tami_base_src/tami_base.hpp"
#include "../src/tami_graph_src/tami_graph.hpp"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>

PYBIND11_MAKE_OPAQUE(std::vector<int>); // alpha_t, epsilon_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiBase::complex_double>); // frequency_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiBase::g_struct>); // g_prod_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiBase::ft_term>); // ft_terms

//PYBIND11_MAKE_OPAQUE(std::vector<TamiGraph::graph_group>); //gg_vec_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiGraph::gg_vec_t>); //gg_matrix_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiGraph::graph_t>); //for TamiGraph::graph_group::graph_vec


namespace py = pybind11;

void init_pytami_wrapper(py::module &m){

    py::bind_vector<std::vector<int>>(m, "VectorInt"); // alpha_t and epsilon_t
    py::bind_vector<std::vector<TamiBase::complex_double>>(m, "frequency_t");
    py::bind_vector<std::vector<TamiBase::g_struct>>(m, "g_prod_t");
    py::bind_vector<std::vector<TamiBase::ft_term>>(m, "ft_terms");

    //py::bind_vector<std::vector<TamiGraph::graph_group>>(m, "gg_vec_t");

    py::class_<TamiBase> TamiBase(m, "TamiBase");
    TamiBase.def(py::init<>());
    TamiBase.def(py::init<at::Device &>()); // device c'tor //TODO: add a  batch_size c'tor
    TamiBase.def("getDevice", &TamiBase::getDevice, "Returns the device stored in the options object.");

    py::class_<TamiBase::ami_vars> (TamiBase, "ami_vars")
        .def(py::init<>())
        .def(py::init<TamiBase::energy_t, TamiBase::frequency_t>())
        .def(py::init<TamiBase::energy_t, TamiBase::frequency_t, double>())
        .def(py::init<TamiBase::energy_t, TamiBase::frequency_t, double, double>())
        .def_readwrite("energy_", &TamiBase::ami_vars::energy_)
        .def_readwrite("frequency_", &TamiBase::ami_vars::frequency_)
        .def_readwrite("prefactor", &TamiBase::ami_vars::prefactor)
        .def_readwrite("BETA_", &TamiBase::ami_vars::BETA_)
        .def_readwrite("gamma_", &TamiBase::ami_vars::gamma_);

    py::class_<TamiBase::ami_parms> (TamiBase, "ami_parms")
        .def(py::init<>())
        .def(py::init<int, double>())
        .def(py::init<int, double, TamiBase::graph_type>())
        .def(py::init<int, double, TamiBase::graph_type, TamiBase::int_type, TamiBase::disp_type>())
        .def_readwrite("N_INT_", &TamiBase::ami_parms::N_INT_)
        .def_readwrite("N_EXT_", &TamiBase::ami_parms::N_EXT_)
        .def_readwrite("E_REG_", &TamiBase::ami_parms::E_REG_)
        .def_readwrite("tol_", &TamiBase::ami_parms::tol_)
        .def_readwrite("TYPE_", &TamiBase::ami_parms::TYPE_)
        .def_readwrite("int_type_", &TamiBase::ami_parms::int_type_)
        .def_readwrite("dispersion_", &TamiBase::ami_parms::dispersion_);

    py::class_<TamiBase::g_struct> (TamiBase, "g_struct")
        .def(py::init<TamiBase::epsilon_t, TamiBase::alpha_t, TamiBase::stat_type>())
        .def(py::init<TamiBase::epsilon_t, TamiBase::alpha_t>())
        .def(py::init<>())
        .def_readwrite("eps_", &TamiBase::g_struct::eps_)
        .def_readwrite("alpha_", &TamiBase::g_struct::alpha_)
        .def_readwrite("stat_", &TamiBase::g_struct::stat_)
        .def_readwrite("species_", &TamiBase::g_struct::species_)
        .def_readwrite("eff_stat_", &TamiBase::g_struct::eff_stat_)
        .def_readwrite("pp", &TamiBase::g_struct::pp);

    py::enum_<TamiBase::graph_type>(TamiBase, "graph_type")
        .value("Sigma", TamiBase::graph_type::Sigma)
        .value("Pi_phuu", TamiBase::graph_type::Pi_phuu)
        .value("Pi_phud", TamiBase::graph_type::Pi_phud)
        .value("Hartree", TamiBase::graph_type::Hartree)
        .value("Bare", TamiBase::graph_type::Bare)
        .value("Greens", TamiBase::graph_type::Greens)
        .value("density", TamiBase::graph_type::density)
        .value("doubleocc", TamiBase::graph_type::doubleocc)
        .value("Pi_ppuu", TamiBase::graph_type::Pi_ppuu)
        .value("Pi_ppud", TamiBase::graph_type::Pi_ppud)
        .value("DOS", TamiBase::graph_type::DOS)
        .value("ENERGY", TamiBase::graph_type::ENERGY)
        .value("FORCE", TamiBase::graph_type::FORCE)
        .export_values();


    py::class_<TamiBase::ft_term> (TamiBase, "ft_term")
        .def(py::init<>())
        .def(py::init<double, TamiBase::FermiTree::fermi_tree_t, TamiBase::g_prod_t>())
        .def_readwrite("sign_", &TamiBase::ft_term::sign_)
        .def_readwrite("ft_", &TamiBase::ft_term::ft_)
        .def_readwrite("g_prod_", &TamiBase::ft_term::g_prod_);

    TamiBase.def("construct", &TamiBase::construct, "Construction function for fermi-tree construction");
    TamiBase.def("evaluate", &TamiBase::evaluate, "Evaluate function for fermi-tree construction");
    TamiBase.def("pretty_print_ft_terms", &TamiBase::pretty_print_ft_terms, "Prints latex formula for the mathematical expression stored in the ft_terms object provided");
    
    
    // Graph bindings 
    // TODO: should this be a separate submodule?

    py::class_<TamiGraph> TamiGraph(m, "TamiGraph");
    TamiGraph.def(py::init<>());
    TamiGraph.def(py::init<TamiBase::graph_type, int>());
    TamiGraph.def(py::init<TamiBase::graph_type, int, int>());

    py::class_<TamiGraph::graph_group> (TamiGraph, "graph_group")
        .def(py::init<>())
        //.def(py::init<double, TamiBase::FermiTree::fermi_tree_t, TamiBase::g_prod_t>())
        .def_readwrite("order_shift", &TamiGraph::graph_group::order_shift)
        .def_readwrite("graph_vec", &TamiGraph::graph_group::graph_vec);

    py::class_<TamiGraph::trojan_graph>(TamiGraph, "trojan_graph")
        //.def(py::init<>())
        .def(py::init<std::vector<TamiGraph::graph_t>, int>())
        .def_readwrite("dummy_var", &TamiGraph::trojan_graph::dummy_var);

    TamiGraph.def("read_ggmp", py::overload_cast<std::string, TamiGraph::gg_matrix_t &, int>(&TamiGraph::read_ggmp), "Reads graph files into ggm from directory provided upto max_ord order.");
    TamiGraph.def("print_ggm", &TamiGraph::print_ggm, "Prints what is contained in the ggm object provided.");
    TamiGraph.def("ggm_label", &TamiGraph::ggm_label, "labels all the graphs contained in the ggm object from order min and up.");
    TamiGraph.def("graph_to_R0", &TamiGraph::graph_to_R0, "Converts provided graph_t into a TamiBase::g_prod_t object.");
    TamiGraph.def("generate_sigma_ct", &TamiGraph::generate_sigma_ct, "generates all Counter term diagrams for the graph provided and store them in the vector provided upto max_number of insertions.");
    TamiGraph.def("get_prefactor", &TamiGraph::get_prefactor, "returns the multiplicative factor of (-1)^fermionic loops.");

    TamiGraph.def("trojan_graph_to_R0", &TamiGraph::trojan_graph_to_R0, "Workaround to get R0 object from graph object inside trojan_graph onject provided.");
    TamiGraph.def("trojan_get_prefactor", &TamiGraph::trojan_get_prefactor, "Workaround to get prefactor from graph object inside trojan_graph onject provided.");

    py::bind_vector<std::vector<TamiGraph::gg_vec_t>>(TamiGraph, "gg_matrix_t");
    py::bind_vector<std::vector<TamiGraph::graph_t>>(TamiGraph, "graph_vector");
}