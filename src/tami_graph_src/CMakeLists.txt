# Library: tami_graph
add_library(tami_graph)

target_sources(tami_graph PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/tamigraph_helpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tamigraph_labelling.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tamigraph_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tamigraph_constructs.cpp
)

target_include_directories(tami_graph
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Graph depends on base; propagate includes transitively
target_link_libraries(tami_graph
  PUBLIC
    tami_base      # brings Boost include dirs via tami_base
    Boost::boost   # (keep PUBLIC if graph headers use Boost directly)
  PRIVATE
    Torch::Torch
)


install(TARGETS tami_graph
  EXPORT torchamiTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/tami_graph_src
  FILES_MATCHING PATTERN "*.hpp"
)
