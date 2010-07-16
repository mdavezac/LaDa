set(PYVASPEXTRACT_PY _dft.py  _gw.py  __init__.py)
install(FILES ${PYVASPEXTRACT_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vasp/extract)

