set( PYCRYSTAL_PY  gruber.py  __init__.py A2BX4.py bravais.py binary.py 
     read_write.py _docstring.py ABX3.py _modify_atom_classes.py)

set(PYCRYSTAL_SOURCES module.cc)

set(PYATOM_HEADERS atom/cdi.hpp atom/type.hpp atom/getset.hpp atom/sequence.hpp
                   atom/sequence_iterator.hpp atom/members.h atom/atom.h)
set(PYCRYSTAL_SOURCES module.cc)


include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})

add_library (pycrystal SHARED ${PYATOM_HEADERS} ${PYCRYSTAL_SOURCES})
target_link_libraries(pycrystal crystal ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY})
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pycrystal PROPERTIES PREFIX "")
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)
# configure_file (
#  "${PROJECT_SOURCE_DIR}/crystal/python/wrappers_init.py"
#  "${PROJECT_BINARY_DIR}/crystal/python/__init__.py"
# )
install( FILES ${PROJECT_BINARY_DIR}/crystal/python/__init__.py DESTINATION
         ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)
if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
