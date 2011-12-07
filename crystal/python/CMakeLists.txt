set( PYCRYSTAL_PY  gruber.py  __init__.py A2BX4.py bravais.py binary.py 
     read_write.py _docstring.py ABX3.py _modify_atom_classes.py)
set( PYCRYSTAL_HEADERS 
     structure.hpp add_atom.hpp structure_contains.hpp atom/atom.h )
   #    atom.hpp  enumerate.hpp  lattice.hpp  layerdepth.hpp  neighbors.hpp
   #    read_structure.hpp  smith.hpp  structure.hpp  symmetry_operator.hpp
   #    which_site.hpp misc.hpp periodic_dnc.hpp )

set( PYCRYSTAL_SOURCES
     python.cc )# structure.cc )
   #    atom.cc  enumerate.cc  lattice.cc  layerdepth.cc  neighbors.cc  python.cc
   #    read_structure.cc  smith.cc  structure.cc  symmetry_operator.cc
   #    which_site.cc misc.cc periodic_dnc.cc )

set(PYATOM_HEADERS atom/cdi.hpp atom/type.hpp atom/getset.hpp atom/set.hpp
                   atom/set_iterator.hpp atom/extract_species.h )
set(PYATOM_SOURCES atom/module.cc)


include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
python_add_module (pycrystal ${PYCRYSTAL_HEADERS} ${PYCRYSTAL_SOURCES})
target_link_libraries(pycrystal crystal ${Boost_PYTHON_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pycrystal)
mark_as_advanced(PYTHON_MODULE_pycrystal_BUILD_SHARED)
set_target_properties(pycrystal PROPERTIES PREFIX "_")
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME crystal)
install(FILES ${PYCRYSTAL_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

add_library (pyatom SHARED ${PYATOM_HEADERS} ${PYATOM_SOURCES})
target_link_libraries(pyatom crystal ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} )
set_target_properties(pyatom PROPERTIES OUTPUT_NAME atom)
set_target_properties(pyatom PROPERTIES PREFIX "")
install(TARGETS pyatom DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal/cppwrappers)
configure_file (
 "${PROJECT_SOURCE_DIR}/crystal/python/wrappers_init.py"
 "${PROJECT_BINARY_DIR}/crystal/python/__init__.py"
)
install( FILES ${PROJECT_BINARY_DIR}/crystal/python/__init__.py DESTINATION
         ${CMAKE_PYINSTALL_PREFIX}/lada/crystal/cppwrappers)
# python_add_module (pyatom ${PYATOM_HEADERS} ${PYATOM_SOURCES})
# target_link_libraries(pyatom crystal ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} )
# mark_as_advanced(PYTHON_ENABLE_MODULE_pyatom)
# mark_as_advanced(PYTHON_MODULE_pyatom_BUILD_SHARED)
# set_target_properties(pyatom PROPERTIES OUTPUT_NAME atom)
# install(TARGETS pyatom DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal/cppwrappers)
# configure_file (
#  "${PROJECT_SOURCE_DIR}/crystal/python/wrappers_init.py"
#  "${PROJECT_BINARY_DIR}/crystal/python/__init__.py"
# )
# install( FILES ${PROJECT_BINARY_DIR}/crystal/python/__init__.py DESTINATION
#          ${CMAKE_PYINSTALL_PREFIX}/lada/crystal/cppwrappers)

if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
