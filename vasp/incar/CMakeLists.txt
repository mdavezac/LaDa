set(PYVASPINCAR_PY __init__.py   _params.py)
install(FILES ${PYVASPINCAR_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/vasp/incar
                PERMISSIONS ${PERMISSIONS})

