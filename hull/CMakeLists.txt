set(PYHULL_PY hull2d.py hull3d.py __init__.py)

install(FILES ${PYHULL_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/hull
                PERMISSIONS ${PERMISSIONS})
