set(PYCE_PY _fit.py __init__.py)
set( PYCE_HEADERS 
     ce.hpp  ce.impl.hpp  clusters.hpp  create_clusters.hpp create_pairs.hpp
     find_pis.hpp  mlcluster.hpp  mlclusters.hpp )
set( PYCE_SOURCES
     ce.cc  clusters.cc  create_clusters.cc  create_pairs.cc  find_pis.cc
     mlcluster.cc  mlclusters.cc  python.cc )


include_directories(${PYTHON_INCLUDE_DIRS})
python_add_module (pyce ${PYCE_HEADERS} ${PYCE_SOURCES})
target_link_libraries(pyce ce ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyce)
mark_as_advanced(PYTHON_MODULE_pyopt_BUILD_SHARED)
set_target_properties(pyce PROPERTIES PREFIX "_")
set_target_properties(pyce PROPERTIES OUTPUT_NAME ce)
install(FILES ${PYCE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ce)
install(TARGETS pyce DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ce)
