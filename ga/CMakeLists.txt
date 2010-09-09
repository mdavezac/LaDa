set(PYGA_PY darwin.py random_search.py  standard.py bitstring.py)

set(which_gapackages "'bitstring', 'standard', 'darwin', 'random_search'")

if(ce_functional)
  set(PYGA_PY ${PYGA_PY} ce.py)
  set(which_gapackages "${which_gapackages}, 'ce'")
endif(ce_functional)


if(escan_functional)
  add_subdirectory(escan)
  set(which_gapackages "${which_gapackages}, 'escan'")
endif(escan_functional)


configure_file (
 "${PROJECT_SOURCE_DIR}/ga/__init__.in.py"
 "${PROJECT_BINARY_DIR}/ga/__init__.py"
)
install(FILES ${PYGA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga)
install(FILES ${PROJECT_BINARY_DIR}/ga/__init__.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga)
