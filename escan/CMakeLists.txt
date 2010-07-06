set( PYESCAN_PY
     _bandgap.py  band_structure.py  derivatives.py  _extract.py  __init__.py
     _methods.py  va.py  _wfns.py )
set(PYESCAN_F90 escan_types.f90 fortran.f90 wavefunctions_module.f90)
set(PYESCAN_HEADERS escan.hpp  wavefunctions.hpp)
set(PYESCAN_SOURCES escan.cc python.cc wavefunctions.cc)
include_directories(MPI_INCLUDE_PATH)

python_add_module (pyescan ${PYESCAN_F90} ${PYESCAN_HEADERS} ${PYESCAN_SOURCES})

include_directories(${MPI_INCLUDE_PATH} ${PYTHON_INCLUDE_DIRS} ${ESCAN_INCLUDE_DIRS})
target_link_libraries(pyescan ${Boost_PYTHON_LIBRARY} ${MPI_LIBRARIES} ${ESCAN_LIBRARY} ${GENPOT_LIBARY})
mark_as_advanced(PYTHON_ENABLE_MODULE_pyescan)
mark_as_advanced(PYTHON_MODULE_pyescan_BUILD_SHARED)
set_target_properties(pyescan PROPERTIES PREFIX "_")
set_target_properties(pyescan PROPERTIES OUTPUT_NAME escan)

install(TARGETS pyescan DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/escan)
install(FILES ${PYESCAN_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/escan)
