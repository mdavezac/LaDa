set(PYPT_PY __init__.py _create_data.py _element.py)

install(FILES ${PYPT_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/periodic_table)
