#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_pytami_wrapper(py::module &);

namespace pytami{

    PYBIND11_MODULE(pytami, m){
        // optional docstring
        m.doc() = "Pytami: a torch enabled python library to evaluate temporal integrals which arise in Feynman diagrams by Algorithmic Matsubara Integration."; // TODO: maybe change this?
        init_pytami_wrapper(m);
    }
}