set(PYCE_PY _fit.py)
set( PYCE_HEADERS 
     ce.hpp  ce.impl.hpp  clusters.hpp  create_clusters.hpp create_pairs.hpp
     find_pis.hpp  mlcluster.hpp  mlclusters.hpp )
set( PYCE_SOURCES
     ce.cc  clusters.cc  create_clusters.cc  create_pairs.cc  find_pis.cc
     mlcluster.cc  mlclusters.cc  python.cc )


include_directories(${PYTHON_INCLUDE_DIRS} ${PYUBLAS_INCLUDE_DIRS})
python_add_module (pyce ${PYCE_HEADERS} ${PYCE_SOURCES})
target_link_libraries(pyce ce ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyce)
mark_as_advanced(PYTHON_MODULE_pyce_BUILD_SHARED)
set_target_properties(pyce PROPERTIES PREFIX "_")
set_target_properties(pyce PROPERTIES OUTPUT_NAME ce)

install(FILES ${PYCE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ce)
install(TARGETS pyce DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ce)

set(LATTICE_TYPE "cubic")
if(${LATTICE_TYPE} STREQUAL "cubic")
  add_definitions(-D_CUBIC_CE_)
elseif(${LATTICE_TYPE} STREQUAL "tetra")
  add_definitions(-D_TETRA_CE_)
endif(${LATTICE_TYPE} STREQUAL "cubic")
configure_file (
 "${PROJECT_SOURCE_DIR}/ce/python/__init__.in.py"
 "${PROJECT_BINARY_DIR}/ce/python/__init__.py"
)
install(FILES ${PROJECT_BINARY_DIR}/ce/python/__init__.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ce)

