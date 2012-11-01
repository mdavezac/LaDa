namespace LaDa
{
  namespace math
  {
    static PyObject* is_integer(PyObject *_module, PyObject* _in)
    {
      using namespace LaDa;
      if(not PyArray_Check(_in))
      {
        LADA_PYERROR(TypeError, "Argument should be a numpy array.");
        return NULL;
      }
      int const type = PyArray_TYPE(_in);
#     ifdef  LADA_NPYITER
#       error LADA_NPYITER already defined
#     endif
#     define LADA_NPYITER(TYPE, NUM_TYPE)                                     \
        if(type == NUM_TYPE)                                                  \
        {                                                                     \
          python::Object iterator = PyArray_IterNew(_in);                     \
          while(PyArray_ITER_NOTDONE(iterator.borrowed()))                    \
          {                                                                   \
            TYPE const x = *((TYPE*)PyArray_ITER_DATA(iterator.borrowed()));  \
            if(not LaDa::math::eq(x, TYPE(std::floor(x+0.1))) )               \
              { Py_RETURN_FALSE; }                                            \
            PyArray_ITER_NEXT(iterator.borrowed());                           \
          }                                                                   \
          Py_RETURN_TRUE;                                                     \
        }
      LADA_NPYITER( npy_float,      NPY_FLOAT)      
      else LADA_NPYITER( npy_double,     NPY_DOUBLE     )
      else LADA_NPYITER( npy_longdouble, NPY_LONGDOUBLE )
      else if(    type == NPY_INT
               or type == NPY_UINT        
               or type == NPY_LONG        
               or type == NPY_LONGLONG    
               or type == NPY_ULONGLONG   
               or type == NPY_BYTE        
               or type == NPY_SHORT       
               or type == NPY_USHORT     ) Py_RETURN_TRUE;
      else
      {
        LADA_PYERROR(TypeError, "Unknown numpy array type.");
        return NULL;
      }
#     undef LADA_NPYITER
    }

    static PyObject* pyfloor_int(PyObject *_module, PyObject* _in)
    {
      using namespace LaDa;
      if(not PyArray_Check(_in))
      {
        LADA_PYERROR(TypeError, "Argument should be a numpy array.");
        return NULL;
      }
      python::Object result = PyArray_SimpleNew( PyArray_NDIM(_in),
                                                 PyArray_DIMS(_in), 
                                                 NPY_LONG );
      if(not result) return NULL;
      python::Object iter_in = PyArray_IterNew(_in);
      if(not iter_in) return NULL;
      python::Object iter_out = PyArray_IterNew(result.borrowed());
      if(not iter_out) return NULL;
      PyObject* py_iterin = iter_in.borrowed();
      PyObject* py_iterout = iter_out.borrowed();

      int const type = PyArray_TYPE(_in);
#     ifdef  LADA_NPYITER
#       error LADA_NPYITER already defined
#     endif
#     define LADA_NPYITER(TYPE, NUM_TYPE)                                       \
        if(type == NUM_TYPE)                                                    \
        {                                                                       \
          while(PyArray_ITER_NOTDONE(py_iterin))                                \
          {                                                                     \
            *((npy_long*)PyArray_ITER_DATA(py_iterout))                         \
                = math::floor_int<TYPE>(*(TYPE*) PyArray_ITER_DATA(py_iterin)); \
            PyArray_ITER_NEXT(py_iterin);                                       \
            PyArray_ITER_NEXT(py_iterout);                                      \
          }                                                                     \
        }
      LADA_NPYITER( npy_float,      NPY_FLOAT)      
      else LADA_NPYITER( npy_double,     NPY_DOUBLE     )
      else LADA_NPYITER( npy_longdouble, NPY_LONGDOUBLE )
      else if(    type == NPY_INT
               or type == NPY_UINT        
               or type == NPY_LONG        
               or type == NPY_LONGLONG    
               or type == NPY_ULONGLONG   
               or type == NPY_BYTE        
               or type == NPY_SHORT       
               or type == NPY_USHORT     ) Py_RETURN_TRUE;
      else
      {
        LADA_PYERROR(TypeError, "Unknown numpy array type.");
        return NULL;
      }
#     undef LADA_NPYITER
      return result.release();
    }

    static PyObject* Rotation1( PyObject *_module, 
                                PyObject *_args, 
                                PyObject *_kwargs )
    {
      double angle;
      PyObject *_vector;
      static char *kwlist[] = { const_cast<char*>("angle"),
                                const_cast<char*>("direction"), NULL};
      if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "dO:Rotation",
                                          kwlist, &angle, &_vector ) )
          return NULL;
      rVector3d vector;
      if(not python::convert_to_vector(_vector, vector)) return NULL;
      
#     ifndef LADA_WITH_EIGEN3 
        // \typedef type of the affine transformations.
        typedef Eigen::Transform<types::t_real, 3> Affine;
#     else
        // \typedef type of the affine transformations.
        typedef Eigen::Transform<types::t_real, 3, Eigen::Isometry> Affine;
#     endif
      // \typedef type of the angle axis object to initialize roations.
      typedef Eigen::AngleAxis<types::t_real> AngleAxis;
      npy_intp dims[2] = {4, 3};
      int const type = numpy::type<types::t_real>::value;
      PyArrayObject *result = (PyArrayObject*)PyArray_ZEROS(2, dims, type, 1);
      if(not result) return NULL;
      
      rVector3d axis;
      python::convert_to_vector(_vector, axis);
      Affine a;
      a = AngleAxis(angle, axis);
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          *((types::t_real*)(result->data + i*result->strides[0] + j*result->strides[1])) = a(i, j);
      for(size_t j(0); j < 3; ++j)
        *((types::t_real*)(result->data + 3*result->strides[0] + j*result->strides[1])) = 0;
      return (PyObject*)result;
    }
    static PyObject *translation(PyObject *_module, PyObject *_args)
    {
#     ifndef LADA_WITH_EIGEN3 
        // \typedef type of the affine transformations.
        typedef Eigen::Transform<types::t_real, 3> Affine;
#     else
        // \typedef type of the affine transformations.
        typedef Eigen::Transform<types::t_real, 3, Eigen::Isometry> Affine;
#     endif
      // \typedef type of the angle axis object to initialize roations.
      npy_intp dims[2] = {4, 3};
      int const type = numpy::type<types::t_real>::value;
      PyArrayObject *result = (PyArrayObject*)PyArray_ZEROS(2, dims, type, 1);
      if(not result) return NULL;
      
      rVector3d trans;
      python::convert_to_vector(_args, trans);
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          *((types::t_real*)(result->data + i*result->strides[0] + j*result->strides[1])) = i == j? 1: 0;
      for(size_t j(0); j < 3; ++j)
        *((types::t_real*)(result->data + 3*result->strides[0] + j*result->strides[1])) = trans(j);
      return (PyObject*)result;
    }


       
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"is_integer",  is_integer, METH_O,
         "True if the input vector or matrix is integer.\n\n"
         "Takes a numpy array as input.\n" }, 
        {"Rotation",  (PyCFunction)Rotation1, METH_KEYWORDS,
         "Rotation of given angle and axis as a 4x3 symmetry operation.\n\n"
         ":param float angle:\n"
         "   Angle of rotation around the input axis.\n"
         ":param axis:\n"
         "   Axis of rotation.\n"
         ":type axis:\n"
         "   sequence of three numbers.\n" }, 
        {"Translation",  (PyCFunction)translation, METH_O,
         "Translation as a 4x3 symmetry operation.\n\n"
         "Takes a sequence of three numbers as input.\n" }, 
        {"floor_int", (PyCFunction)pyfloor_int, METH_O, 
         "Floors an input array.\n\n"
         "Takes a numpy array as input. Returns an integer array.\n" },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    }; // end of static method table.
  }
}

