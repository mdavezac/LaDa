set(PYJOBS_PY __init__.py  runme.py  templates.py runone.py bleeder.py extract.py
              manipulator.py forwarding_dict.py default_params.py jobdict.py)
install(FILES ${PYJOBS_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada/jobs)

