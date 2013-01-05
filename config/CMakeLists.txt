set(CONFIG_FILES escan.py general.py  ladabase.py  mpi.py  vasp.py ipython.py
    dftcrystal.py deepcopy_rig.py patch_ordered_dict.py)
set(EXTRAS_FILES carver_mpi.py redmesa_mpi.py hopper_mpi.py cx1.py)

install(FILES ${CONFIG_FILES} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada/config)
