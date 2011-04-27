set( PYESCAN_PY
     _bandgap.py  _bandstructure.py  derivatives.py  _extract.py  __init__.py
     _methods.py  va.py  _wfns.py functional.py kescan.py _potential.py kpoints.py
     _massextract.py ldos.py )
set(PYESCAN_HEADERS escan.hpp  wavefunctions.hpp)
set(PYESCAN_SOURCES escan.cc python.cc wavefunctions.cc)
if(ESCAN_FOUND)
  set(PYESCAN_F90 escan_types.f90 fortran.f90 wavefunctions_module.f90)
else(ESCAN_FOUND)
  set(PYESCAN_F90 "")
endif(ESCAN_FOUND)

python_add_module (pyescan ${PYESCAN_F90} ${PYESCAN_HEADERS} ${PYESCAN_SOURCES})

include_directories( ${MPI_INCLUDE_PATH} 
                     ${PYTHON_INCLUDE_DIRS} 
                     ${ESCAN_INCLUDE_DIRS}
                     ${NUMPY_INCLUDES} )
target_link_libraries( pyescan ${Boost_PYTHON_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
if(ESCAN_FOUND)
  target_link_libraries(pyescan ${ESCAN_LIBRARY} ${GENPOT_LIBRARY})
endif(ESCAN_FOUND)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(pyescan ${MPI_LIBRARY} ${Boost_MPI_LIBRARY})
endif (MPI_FOUND)

mark_as_advanced(PYTHON_ENABLE_MODULE_pyescan)
mark_as_advanced(PYTHON_MODULE_pyescan_BUILD_SHARED)
set_target_properties(pyescan PROPERTIES PREFIX "_")
set_target_properties(pyescan PROPERTIES OUTPUT_NAME escan)

install(TARGETS pyescan DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/escan)
install(FILES ${PYESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/escan)
