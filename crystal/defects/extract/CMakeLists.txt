set(PD_PY  __init__.py _single.py _mass.py)

install(FILES ${PD_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/crystal/defects/extract)
