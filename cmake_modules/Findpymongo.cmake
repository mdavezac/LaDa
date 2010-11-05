if(NOT PYTHONINTERP_FOUND)
  find_package(PythonInterp REQUIRED)
endif(NOT PYTHONINTERP_FOUND)

execute_process (
  COMMAND ${PYTHON_EXECUTABLE} -c "import pymongo, os; print os.path.dirname(pymongo.__file__)"
  OUTPUT_VARIABLE pymongo_path
  )
if (pymongo_path)
  SET(PYMONGO_FOUND TRUE)
else (pymongo_path)
  SET(PYMONGO_FOUND FALSE)
endif (pymongo_path)
