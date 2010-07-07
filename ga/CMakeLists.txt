set(PYGA_PY darwin.py  __init__.py  random_search.py  standard.py bitstring.py)
if(ce_functional)
  set(PYGA_PY ${PYGA_PY} ce.py)
endif(ce_functional)
install(FILES ${PYGA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ga)
if(escan_functional)
  add_subdirectory(pescan)
endif(escan_functional)

