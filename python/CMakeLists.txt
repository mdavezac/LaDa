set(PYLADA_HEADERS debug.hpp misc.hpp numpy_types.h std_map.hpp std_map.impl.hpp std_vector.hpp xml.hpp)
set(PYLADA_PY __init__.py)

install(FILES ${PYLADA_PY} DESTINATION ${CMAKE_PYINSTALL_PREFIX}/lada)
if(NOT noinstall_headers)
  install(FILES ${PYLADA_HEADERS} DESTINATION include/LaDa/python)
endif(NOT noinstall_headers)
