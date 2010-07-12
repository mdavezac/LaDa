set(PYJOBS_PY __init__.py  runme.py  templates.py runone.py)
install(FILES ${PYJOBS_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/jobs)

