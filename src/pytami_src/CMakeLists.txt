# Python extension: pytami
# Keep it thin: just the wrapper; link to libs for the implementation.
pybind11_add_module(pytami MODULE
  ${CMAKE_CURRENT_SOURCE_DIR}/pytami_wrapper.cpp
)

target_link_libraries(pytami
  PRIVATE
    tami_graph        # brings in tami_base transitively
    Torch::Torch
)

# Optional: hide symbols by default (pybind11 sets visibility=hidden already on many platforms)
# target_compile_options(pytami PRIVATE -fvisibility=hidden)

# Install to the active site-packages
# Prefer site-arch if defined, else sitelib.
set(_PY_SITE "${Python3_SITEARCH}")
if(NOT _PY_SITE)
  set(_PY_SITE "${Python3_SITELIB}")
endif()

# install(TARGETS pytami
#   LIBRARY DESTINATION "${_PY_SITE}"
# )

# Install into the prefix (e.g., .../install/lib), not system site-packages
install(TARGETS pytami
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
