set(PYVASPNLEP_PY __init__.py mpopt.py systems.py input_generic.py nlep.py system_params.py nlep_defaults.py postprocessing.py plotfit.py plotbs.py)

install(FILES ${PYVASPNLEP_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/vasp/nlep)
