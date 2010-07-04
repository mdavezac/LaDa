set(PYOPT_PY changedir.py  decorators.py  __init__.py  secant.py   tempdir.py)
set(PYOPT_HEADERS convexhull.impl.hpp  errortuple.hpp  redirect.hpp)
set(PYOPT_SOURCES errortuple.cc  python.cc  redirect.cc)

include_directories(${PYTHON_INCLUDE_DIRS})
python_add_module (pyopt ${PYOPT_HEADERS} ${PYOPT_SOURCES})
