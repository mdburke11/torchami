#include "../src/tami_base.hpp"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>

PYBIND11_MAKE_OPAQUE(std::vector<int>); // alpha_t, epsilon_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiBase::complex_double>); // frequency_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiBase::g_struct>); // g_prod_t
PYBIND11_MAKE_OPAQUE(std::vector<TamiBase::ft_term>); // ft_terms


namespace py = pybind11;

void init_pytami_wrapper(py::module &m){

    py::bind_vector<std::vector<int>>(m, "VectorInt"); // alpha_t and epsilon_t
    py::bind_vector<std::vector<TamiBase::complex_double>>(m, "frequency_t");
    py::bind_vector<std::vector<TamiBase::g_struct>>(m, "g_prod_t");
    py::bind_vector<std::vector<TamiBase::ft_term>>(m, "ft_terms");

    py::class_<TamiBase> TamiBase(m, "TamiBase");
    TamiBase.def(py::init<>());
    TamiBase.def(py::init<at::Device &>()); // device c'tor //TODO: add a  batch_size c'tor
    TamiBase.def_readwrite("options", &TamiBase::options); // TODO: remove since a python analog DNE
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

    py::class_<TamiBase::ft_term> (TamiBase, "ft_term")
        .def(py::init<>())
        .def(py::init<double, TamiBase::FermiTree::fermi_tree_t, TamiBase::g_prod_t>())
        .def_readwrite("sign_", &TamiBase::ft_term::sign_)
        .def_readwrite("ft_", &TamiBase::ft_term::ft_)
        .def_readwrite("g_prod_", &TamiBase::ft_term::g_prod_);

    TamiBase.def("construct", &TamiBase::construct, "Construction function for fermi-tree construction");
    TamiBase.def("evaluate", &TamiBase::evaluate, "Evaluate function for fermi-tree construction");
    TamiBase.def("pretty_print_ft_terms", &TamiBase::pretty_print_ft_terms, "Prints latex formula for the mathematical expression stored in the ft_terms object provided");
    

    

    }