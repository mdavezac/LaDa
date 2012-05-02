set(PYIPY_PY __init__.py _goto.py _explore.py _showme.py _doc.py export.py record.py qstat.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ipython
                PERMISSIONS ${PERMISSIONS})
add_subdirectory(launch)

