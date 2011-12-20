set(CONFIG_FILES escan.py  ipython.py  ladabase.py  mpi.py  vasp.py)
set(EXTRAS_FILES carver_mpi.py redmesa_mpi.py hopper_mpi.py)

install(FILES ${CONFIG_FILES} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/config)
