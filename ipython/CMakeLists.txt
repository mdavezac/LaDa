set(PYIPY_PY __init__.py explore.py savejobs.py goto.py listjobs.py
             showme.py ) # export.py qstat.py)
install(FILES ${PYIPY_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/ipython)
add_subdirectory(launch)

