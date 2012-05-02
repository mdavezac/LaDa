set(LADABASE_PY __init__.py ipython.py vasp.py misc.py filter.py extracted.py
                fere.py massextract.py mu_data.py _local_tar.py push.py)
install(FILES ${LADABASE_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ladabase
                PERMISSIONS ${PERMISSIONS})
