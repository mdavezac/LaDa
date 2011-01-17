set(PYCRYSTAL_PY  gruber.py  __init__.py A2BX4.py bravais.py binary.py read_write.py _docstring.py)
set( PYCRYSTAL_HEADERS 
     atom.hpp  enumerate.hpp  lattice.hpp  layerdepth.hpp  neighbors.hpp
     read_structure.hpp  smith.hpp  structure.hpp  symmetry_operator.hpp
     which_site.hpp )

set( PYCRYSTAL_SOURCES
     atom.cc  enumerate.cc  lattice.cc  layerdepth.cc  neighbors.cc  python.cc
     read_structure.cc  smith.cc  structure.cc  symmetry_operator.cc
     which_site.cc )


include_directories(${PYTHON_INCLUDE_DIRS})
python_add_module (pycrystal ${PYCRYSTAL_HEADERS} ${PYCRYSTAL_SOURCES})
target_link_libraries(pycrystal crystal ${Boost_PYTHON_LIBRARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pycrystal)
mark_as_advanced(PYTHON_MODULE_pycrystal_BUILD_SHARED)
set_target_properties(pycrystal PROPERTIES PREFIX "_")
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME crystal)
install(FILES ${PYCRYSTAL_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
