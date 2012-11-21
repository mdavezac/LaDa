namespace LaDa
{
  namespace python
  {
#   define LADA_CRYSTAL_MODULE 0
#   include "../crystal/python/numpy_types.h"
#   include "../crystal/python/wrap_numpy.h"
#   undef LADA_CRYSTAL_MODULE 
  }
  namespace math
  {
    static PyObject* pyis_integer(PyObject *_module, PyObject* _in)
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
          PyObject* iterator = PyArray_IterNew(_in);                          \
          if(not iterator) return NULL;                                       \
          while(PyArray_ITER_NOTDONE(iterator))                               \
          {                                                                   \
            TYPE const x = *((TYPE*)PyArray_ITER_DATA(iterator));             \
            if(not LaDa::math::eq(x, TYPE(std::floor(x+0.1))) )               \
              { Py_DECREF(iterator); Py_RETURN_FALSE; }                       \
            PyArray_ITER_NEXT(iterator);                                      \
          }                                                                   \
          Py_DECREF(iterator);                                                \
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
      int const type = PyArray_TYPE(_in);
      if(    type == NPY_INT
          or type == NPY_UINT        
          or type == NPY_LONG        
          or type == NPY_LONGLONG    
          or type == NPY_ULONGLONG   
          or type == NPY_BYTE        
          or type == NPY_SHORT       
          or type == NPY_USHORT     )
      {
        LADA_PYERROR(TypeError, "Numpy array is already an interger type.");
        Py_RETURN_TRUE; 
      }

      PyObject* result = PyArray_SimpleNew( PyArray_NDIM(_in),
                                            PyArray_DIMS(_in), 
                                            NPY_LONG );
      if(not result) return NULL;
      PyObject* iter_in = PyArray_IterNew(_in);
      if(not iter_in) {Py_DECREF(result); return NULL;}
      PyObject* iter_out = PyArray_IterNew(result);
      if(not iter_out) {Py_DECREF(iter_in); Py_DECREF(result); return NULL;}

#     ifdef  LADA_NPYITER
#       error LADA_NPYITER already defined
#     endif
#     define LADA_NPYITER(TYPE, NUM_TYPE)                                       \
        if(type == NUM_TYPE)                                                    \
        {                                                                       \
          while(PyArray_ITER_NOTDONE(iter_in))                                  \
          {                                                                     \
            *((npy_long*)PyArray_ITER_DATA(iter_out))                           \
                = math::floor_int<TYPE>(*(TYPE*) PyArray_ITER_DATA(iter_in));   \
            PyArray_ITER_NEXT(iter_in);                                         \
            PyArray_ITER_NEXT(iter_out);                                        \
          }                                                                     \
          Py_DECREF(iter_in);                                                   \
          Py_DECREF(iter_out);                                                  \
        }
      LADA_NPYITER( npy_float,      NPY_FLOAT)      
      else LADA_NPYITER( npy_double,     NPY_DOUBLE     )
      else LADA_NPYITER( npy_longdouble, NPY_LONGDOUBLE )
      else
      {
        Py_DECREF(result);
        LADA_PYERROR(TypeError, "Unknown numpy array type.");
        return NULL;
      }
#     undef LADA_NPYITER
      return result;
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
      int const type = python::numpy::type<types::t_real>::value;
      PyArrayObject *result = (PyArrayObject*)PyArray_ZEROS(2, dims, type, 1);
      if(not result) return NULL;
      
      Affine a;
      a = AngleAxis(angle, vector);
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
      int const type = python::numpy::type<types::t_real>::value;
      PyArrayObject *result = (PyArrayObject*)PyArray_ZEROS(2, dims, type, 1);
      if(not result) return NULL;
      
      rVector3d trans;
      if(not python::convert_to_vector(_args, trans)) return NULL;
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          *((types::t_real*)(result->data + i*result->strides[0] + j*result->strides[1])) = i == j? 1: 0;
      for(size_t j(0); j < 3; ++j)
        *((types::t_real*)(result->data + 3*result->strides[0] + j*result->strides[1])) = trans(j);
      return (PyObject*)result;
    }

    static PyObject* pygruber(PyObject* _module, PyObject* _args, PyObject *_kwargs)
    {
      PyObject *_cell;
      int itermax = 0;
      double tolerance = types::tolerance;
      static char *kwlist[] = { const_cast<char*>("cell"),
                                const_cast<char*>("itermax"),
                                const_cast<char*>("tolerance"), NULL};
      if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "O|id:gruber",
                                          kwlist, &_cell, &itermax, &tolerance ) )
          return NULL;
      rMatrix3d cell;
      if(not python::convert_to_matrix(_cell, cell)) return NULL;

      try
      {
        rMatrix3d result = gruber(cell, itermax, tolerance);
        return python::wrap_to_numpy(result);
      }
      catch(error::singular_matrix& _e)
      {
        LADA_PYERROR(input, "Singular matrix in gruber.");
        return NULL;
      }
      catch(error::infinite_loop& _e)
      {
        LADA_PYERROR(internal, "Maximum number of iterations reached in gruber.");
        return NULL;
      }
      catch(...)
      {
        LADA_PYERROR(internal, "Error encoutered in smith normal form.");
        return NULL;
      }
    }

    static PyObject* pysmith(PyObject* _module, PyObject* _matrix)
    {
      iMatrix3d matrix;
      if(not python::convert_to_matrix(_matrix, matrix)) return NULL;

      iMatrix3d S, L, R;
      try { smith_normal_form(S, L, matrix, R); }
      catch(error::singular_matrix& _e)
      {
        LADA_PYERROR(input, "Cannot compute smith normal form of singular matrix.");
        return NULL;
      }
      catch(...)
      {
        LADA_PYERROR(internal, "Error encoutered in smith normal form.");
        return NULL;
      }
      PyObject *result = PyTuple_New(3);
      if(not result) return NULL;
      PyObject *pyS = python::wrap_to_numpy(S);
      if(not pyS) { Py_DECREF(result); return NULL; }
      PyObject *pyL = python::wrap_to_numpy(L);
      if(not pyL) { Py_DECREF(result); Py_DECREF(pyS); return NULL; }
      PyObject *pyR = python::wrap_to_numpy(R);
      if(not pyR) { Py_DECREF(result); Py_DECREF(pyS); Py_DECREF(pyL); return NULL; }
      PyTuple_SET_ITEM(result, 0, pyS);
      PyTuple_SET_ITEM(result, 1, pyL);
      PyTuple_SET_ITEM(result, 2, pyR);
      return result;
    }
       
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"is_integer",  pyis_integer, METH_O,
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
        {"gruber", (PyCFunction)pygruber, METH_KEYWORDS,
         "Determines Gruber cell of an input cell.\n\n"
         "The Gruber cell is an optimal parameterization of a lattice, eg shortest\n"
         "cell-vectors and angles closest to 90 degrees.\n\n"
         ":param cell:\n  The input lattice cell-vectors.\n"
         ":type cell:  numpy 3x3 array\n"
         ":param int itermax:\n  Maximum number of iterations. Defaults to 0, ie infinite.\n"
         ":param float tolerance:\n  Tolerance parameter when comparing real numbers. "
            "Defaults to LaDa internals.\n"
         ":returns: An equivalent standardized cell.\n"
         ":raises lada.error.input: If the input matrix is singular.\n"
         ":raises lada.error.internal: If the maximum number of iterations is reached.\n"},
         {"smith_normal_form", (PyCFunction)pysmith, METH_O,
          "Computes smith normal form of a matrix.\n\n"
          "If ``M`` is the input matrix, then the smith normal form is\n"
          ":math:`S = L\\cdot M \\cdot R`, with ``S`` diagonal and ``L``\n"
          "and ``R`` are diagonal.\n\nThis function takes an *integer* matrix\n"
          "as its only input.\n\n"
          ":returns: the tuple of matrices (``S``, ``L``, ``R``)." },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    }; // end of static method table.
  }
}

