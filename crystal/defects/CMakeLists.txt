set(PD_PY  __init__.py extract.py)

if(vasp_functional AND jobs_package)
  install(FILES ${PD_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/crystal/defects)
endif(vasp_functional AND jobs_package)
