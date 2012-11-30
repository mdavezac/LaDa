install(FILES error.py DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada)
if(NOT noinstall_headers)
  install(FILES exceptions.h DESTINATION include/LaDa/errors)
endif(NOT noinstall_headers)
