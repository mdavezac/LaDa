set(LIBCRYSTAL_HEADERS atom/atom.h atom/pybase.h
                       structure/structure.h structure/pybase.h
                       hart-forcade/hart-forcade.h hart-forcade/pybase.h
                       equivalent_structures.h supercell.h primitive.h
                       space_group.h map_sites.h periodic_dnc.h neighbors.h
                       coordination_shells.h confsplit.h )
set(LIBCRYSTAL_PYHEADERS atom/cdi.hpp atom/getset.hpp atom/members.hpp
                         structure/cdi.hpp structure/getset.hpp structure/sequence.hpp 
                         structure/members.hpp 
                         hart-forcade/cdi.hpp hart-forcade/get.hpp
                         hart-forcade/members.hpp hart-forcade/macro.hpp )
set(LIBCRYSTAL_SOURCES atom/pybase.cc 
                       structure/pybase.cc structure/structure.cc
                       hart-forcade/pybase.cc
                       utilities.cc primitive.cc equivalent_structures.cc
                       space_group.cc periodic_dnc.cc map_sites.cc
                       neighbors.cc coordination_shells.cc confsplit.cc )
set(LIBCRYSTAL_PYSOURCES __init__.py read.py write.py ABX.py A2BX4.py bravais.py binary.py)

include_directories(${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDES})
add_library (crystal SHARED ${LIBCRYSTAL_HEADERS} ${LIBCRYSTAL_PYHEADERS} ${LIBCRYSTAL_SOURCES})
target_link_libraries(crystal ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARY} mathlib)
install(TARGETS crystal DESTINATION lib64)

add_library (pycrystal SHARED module.cc)
target_link_libraries(pycrystal ${PYTHON_LIBRARY} crystal)
set_target_properties(pycrystal PROPERTIES OUTPUT_NAME cppwrappers)
set_target_properties(pycrystal PROPERTIES PREFIX "")
install(TARGETS pycrystal DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)
install(FILES ${LIBCRYSTAL_PYSOURCES} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal)

if(compile_tests OR compile_debug_tests)
  add_subdirectory(tests)
endif(compile_tests OR compile_debug_tests)
