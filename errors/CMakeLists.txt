install(FILES error.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/pylada)
if(NOT noinstall_headers)
  install(FILES exceptions.h DESTINATION include/Pylada/errors)
endif(NOT noinstall_headers)
