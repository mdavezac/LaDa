namespace LaDa
{
  namespace crystal
  {
#   ifdef LADA_MACRO 
#     error LADA_MACRO ALREADY defined.
#   endif
#   define LADA_MACRO(NAME)                                                             \
    /** Wrapper around NAME function */                                                 \
    PyObject* NAME ## _wrapper(PyObject *_module, PyObject *_args)                      \
    {                                                                                   \
      Py_ssize_t const N = PyTuple_Size(_args);                                         \
      if(N != 3 and N != 2)                                                             \
      {                                                                                 \
        LADA_PYERROR(TypeError, #NAME " expects two vectors and a matrix as input.");   \
        return NULL;                                                                    \
      }                                                                                 \
      math::rMatrix3d cell, invcell;                                                    \
      math::rVector3d a;                                                                \
      if(not python::convert_to_vector(PyTuple_GET_ITEM(_args, 0), a)) return false;    \
      if(not python::convert_to_matrix(PyTuple_GET_ITEM(_args, 1), cell)) return false; \
      if(N == 3 and not python::convert_to_matrix(PyTuple_GET_ITEM(_args, 2), invcell)) \
        return false;                                                                   \
      if(N == 2) invcell = cell.inverse();                                              \
      try                                                                               \
      {                                                                                 \
        math::rVector3d vector = NAME(a, cell, invcell);                                \
        npy_intp dim[1] = {3};                                                          \
        PyObject* result = PyArray_SimpleNew(1, dim,                                    \
            math::numpy::type<math::rVector3d::Scalar>::value );                        \
        if(result == NULL) return NULL;                                                 \
        char * const elem = (char*const) PyArray_DATA(result);                          \
        npy_intp const stride = PyArray_STRIDE(result, 0);                              \
        *((math::rVector3d::Scalar*)elem) = vector[0];                                  \
        *((math::rVector3d::Scalar*)(elem+stride)) = vector[1];                         \
        *((math::rVector3d::Scalar*)(elem+(stride<<1))) = vector[2];                    \
        return result;                                                                  \
      }                                                                                 \
      /* catch exceptions. */                                                           \
      catch(error::pyerror &e)                                                          \
      {                                                                                 \
        /* Returns if python error already set. */                                      \
        if(PyErr_Occurred() != NULL) return NULL;                                       \
        /* Otherwise, throw our own */                                                  \
        std::string const * message = boost::get_error_info<error::string>(e);          \
        PyObject ** exception = boost::get_error_info<error::pyexcept>(e);              \
        if(exception == NULL)                                                           \
        {                                                                               \
          LADA_PYERROR(InternalError, message ? message->c_str():                       \
                           ( "Caught exception in " #NAME ": " +                        \
                             std::string(e.what()) ).c_str() );                         \
        }                                                                               \
        else if(message == NULL)                                                        \
          PyErr_SetString( *exception,                                                  \
                           ( "Caught exception in " #NAME ": " +                        \
                             std::string(e.what()) ).c_str() );                         \
        else PyErr_SetString(*exception, message->c_str());                             \
      }                                                                                 \
      catch(std::exception &e)                                                          \
      {                                                                                 \
        std::string const message =  "Caught exception in " #NAME                       \
                                     + std::string(e.what());                           \
        LADA_PYERROR(InternalError, message.c_str());                                   \
      }                                                                                 \
      catch(...)                                                                        \
      {                                                                                 \
        LADA_PYERROR(InternalError, "Caught unknown c++ exception in " #NAME ".");      \
      }                                                                                 \
      return NULL;                                                                      \
    }
    LADA_MACRO(into_cell)
    LADA_MACRO(into_voronoi)
    LADA_MACRO(zero_centered)
#   undef into_cell
    //! Wrapper around periodic image tester.
    PyObject* are_periodic_images_wrapper(PyObject *_module, PyObject *_args)
    {
      Py_ssize_t const N = PyTuple_Size(_args);
      if(N != 3 and N != 4)
      {
        LADA_PYERROR(TypeError, "are_periodic_images expects two vectors and a matrix as input.");
        return NULL;
      }
      math::rMatrix3d invcell;
      math::rVector3d a, b;
      if(not python::convert_to_vector(PyTuple_GET_ITEM(_args, 0), a)) return false;
      if(not python::convert_to_vector(PyTuple_GET_ITEM(_args, 1), b)) return false;
      if(not python::convert_to_matrix(PyTuple_GET_ITEM(_args, 2), invcell)) return false;
      try
      { 
        if(N == 3 and math::are_periodic_images(a, b, invcell) ) Py_RETURN_TRUE;
        else if(N == 4)
        {
          types::t_real tolerance = types::tolerance;
          PyObject* o = PyTuple_GET_ITEM(_args, 3);
          if(PyFloat_Check(o)) tolerance = PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o)) tolerance = PyInt_AS_LONG(o);
          else
          {
            LADA_PYERROR(TypeError, "Tolerance should be a number.");
            return NULL;
          }
          if(math::are_periodic_images(a, b, invcell, tolerance)) Py_RETURN_TRUE;
        }
        Py_RETURN_FALSE;
      }
      // catch exceptions.
      catch(error::pyerror &e)
      {
        // Returns if python error already set.
        if(PyErr_Occurred() != NULL) return NULL;
        // Otherwise, throw our own
        std::string const * message = boost::get_error_info<error::string>(e);
        PyObject ** exception = boost::get_error_info<error::pyexcept>(e);
        if(exception == NULL)
        {
          LADA_PYERROR(InternalError, message ? message->c_str(): 
                           ( "Caught exception while testing periodic images: " + 
                             std::string(e.what()) ).c_str() );
        }
        else if(message == NULL)
          PyErr_SetString( *exception, 
                           ( "Caught exception while testing periodic images: " + 
                             std::string(e.what()) ).c_str() );
        else PyErr_SetString(*exception, message->c_str());
      }
      catch(std::exception &e)
      {
        std::string const message =  "Caught exception while testing periodic images: "
                                     + std::string(e.what());
        LADA_PYERROR(InternalError, message.c_str()); 
      }
      catch(...)
      { 
        LADA_PYERROR(InternalError, "Caught unknown c++ exception in are_periodic_images.");
      }
      return NULL;
    }
    //! Wrapper around the supercell method.
    PyObject* supercell_wrapper(PyObject *_module, PyObject *_args, PyObject *_kwargs)
    {
      // check/convert input parameters.
      PyObject* lattice;
      PyObject* cellin;
      static char *kwlist[] = { const_cast<char*>("lattice"),
                                const_cast<char*>("supercell"), NULL};
      if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "OO:supercell", kwlist,
                                          &lattice, &cellin) )
        return NULL;
      if(not PyStructure_Check(lattice))
      {
        LADA_PYERROR(TypeError, "Input is not a crystal.Structure object."); 
        return NULL;
      }
      math::rMatrix3d cell;
      if(not python::convert_to_matrix(cellin, cell)) return NULL;
      // create supercell.
      try { return supercell(Structure::acquire(lattice), cell).release(); } 
      // catch exceptions.
      catch(error::pyerror &e)
      {
        // Returns if python error already set.
        if(PyErr_Occurred() != NULL) return NULL;
        // Otherwise, throw our own
        std::string const * message = boost::get_error_info<error::string>(e);
        PyObject ** exception = boost::get_error_info<error::pyexcept>(e);
        if(exception == NULL)
        {
          LADA_PYERROR(InternalError, message ? message->c_str(): 
                           ( "Caught exception while creating supercell: " + 
                             std::string(e.what()) ).c_str() );
        }
        else if(message == NULL)
          PyErr_SetString( *exception, 
                           ( "Caught exception while creating supercell: " + 
                             std::string(e.what()) ).c_str() );
        else PyErr_SetString(*exception, message->c_str());
      }
      catch(std::exception &e)
      {
        std::string const message =  "Caught exception while creating supercell: "
                                     + std::string(e.what());
        LADA_PYERROR(InternalError, message.c_str()); 
      }
      catch(...) { LADA_PYERROR(InternalError, "Caught unknown c++ exception."); }
      return NULL;
    }

    PyObject* primitive_wrapper(PyObject* _module, PyObject *_args)
    {
      Py_ssize_t const N(PyTuple_Size(_args));
      if(N != 1 and N != 2)
      {
        LADA_PYERROR(TypeError, "primitive expects a structure as input, and optionally a tolerance.");
        return NULL;
      }
      if(not PyStructure_Check(PyTuple_GET_ITEM(_args, 0)))
      {
        LADA_PYERROR(TypeError, "First input argument to primitive is not a structure.");
        return NULL;
      }
      Structure const lattice = Structure::acquire(PyTuple_GET_ITEM(_args, 0));
      PyObject *tolobj = N == 2? PyTuple_GET_ITEM(_args, 1): NULL;
      if(N == 2 and (PyInt_Check(tolobj) == false and PyFloat_Check(tolobj) == false))
      {
        LADA_PYERROR(TypeError, "First input argument to primitive is not a structure.");
        return NULL;
      }
      types::t_real tolerance = N == 1 ? -1:
              ( PyInt_Check(tolobj) ? PyInt_AS_LONG(tolobj): PyFloat_AS_DOUBLE(tolobj) );
      try { return primitive(lattice, tolerance).release(); } 
      // catch exceptions.
      catch(error::pyerror &e)
      {
        // Returns if python error already set.
        if(PyErr_Occurred() != NULL) return NULL;
        // Otherwise, throw our own
        std::string const * message = boost::get_error_info<error::string>(e);
        PyObject ** exception = boost::get_error_info<error::pyexcept>(e);
        if(exception == NULL)
        {
          LADA_PYERROR(InternalError, message ? message->c_str(): 
                           ( "Caught exception while creating primitive lattice: " + 
                             std::string(e.what()) ).c_str() );
        }
        else if(message == NULL)
          PyErr_SetString( *exception, 
                           ( "Caught exception while creating primitive structure: " + 
                             std::string(e.what()) ).c_str() );
        else PyErr_SetString(*exception, message->c_str());
      }
      catch(std::exception &e)
      {
        std::string const message =  "Caught exception while creating primitive structure: "
                                     + std::string(e.what());
        LADA_PYERROR(InternalError, message.c_str()); 
      }
      catch(...) { LADA_PYERROR(InternalError, "Caught unknown c++ exception."); }
      return NULL;
    };

    PyObject* is_primitive_wrapper(PyObject* _module, PyObject *_args)
    {
      Py_ssize_t const N(PyTuple_Size(_args));
      if(N != 1 and N != 2)
      {
        LADA_PYERROR(TypeError, "primitive expects a structure as input, and optionally a tolerance.");
        return NULL;
      }
      if(not PyStructure_Check(PyTuple_GET_ITEM(_args, 0)))
      {
        LADA_PYERROR(TypeError, "First input argument to primitive is not a structure.");
        return NULL;
      }
      Structure const lattice = Structure::acquire(PyTuple_GET_ITEM(_args, 0));
      PyObject *tolobj = N == 2? PyTuple_GET_ITEM(_args, 1): NULL;
      if(N == 2 and (PyInt_Check(tolobj) == false and PyFloat_Check(tolobj) == false))
      {
        LADA_PYERROR(TypeError, "First input argument to primitive is not a structure.");
        return NULL;
      }
      types::t_real tolerance = N == 1 ? -1:
              ( PyInt_Check(tolobj) ? PyInt_AS_LONG(tolobj): PyFloat_AS_DOUBLE(tolobj) );
      try
      {
        if(is_primitive(lattice, tolerance)) Py_RETURN_TRUE;
        Py_RETURN_FALSE;
      } 
      // catch exceptions.
      catch(error::pyerror &e)
      {
        // Returns if python error already set.
        if(PyErr_Occurred() != NULL) return NULL;
        // Otherwise, throw our own
        std::string const * message = boost::get_error_info<error::string>(e);
        PyObject ** exception = boost::get_error_info<error::pyexcept>(e);
        if(exception == NULL)
        {
          LADA_PYERROR(InternalError, message ? message->c_str(): 
                           ( "Caught exception while testing whether structure is primitive: " + 
                             std::string(e.what()) ).c_str() );
        }
        else if(message == NULL)
          PyErr_SetString( *exception, 
                           ( "Caught exception while testing whether structure is primitive: " + 
                             std::string(e.what()) ).c_str() );
        else PyErr_SetString(*exception, message->c_str());
      }
      catch(std::exception &e)
      {
          std::string const message =  "Caught exception while testing whether structure is primitive: "
                                       + std::string(e.what());
          LADA_PYERROR(InternalError, message.c_str()); 
        }
        catch(...) { LADA_PYERROR(InternalError, "Caught unknown c++ exception."); }
        return NULL;
      };

      PyObject* cell_invariants_wrapper(PyObject *_module, PyObject *_args)
      {
        Py_ssize_t const N(PyTuple_Size(_args));
        if(N != 1 and N != 2)
        {
          LADA_PYERROR(TypeError, "cell_invariants expects one or two arguments only.");
          return NULL;
        }
        PyObject * const arg0 = PyTuple_GET_ITEM(_args, 0);
        math::rMatrix3d cell;
        if(PyStructure_Check(arg0)) cell = ((StructureData*)arg0)->cell;
        else if(not python::convert_to_matrix(arg0, cell)) return NULL;
        types::t_real tolerance = types::tolerance;
        if(N == 2)
        {
          PyObject * const arg1 = PyTuple_GET_ITEM(_args, 1);
          if(PyInt_Check(arg1)) tolerance = PyInt_AS_LONG(arg1);
          else if(PyFloat_Check(arg1)) tolerance = PyFloat_AS_DOUBLE(arg1);
          else
          {
            LADA_PYERROR(TypeError, "Second argument to cell_invariants should a real number.");
            return NULL;
          }
        }
        return cell_invariants(cell, tolerance);
      }
      
      PyObject* space_group_wrapper(PyObject *_module, PyObject *_args)
      {
        Py_ssize_t const N(PyTuple_Size(_args));
        if(N != 1 and N != 2)
        {
          LADA_PYERROR(TypeError, "space_group expects one or two arguments only.");
          return NULL;
        }
        PyObject * const arg0 = PyTuple_GET_ITEM(_args, 0);
        if(not PyStructure_Check(arg0)) 
        {
          LADA_PYERROR(TypeError, "space_group expects a Structure as first argument.");
          return NULL;
        }
        Structure structure = Structure::acquire(arg0); 
        types::t_real tolerance = types::tolerance;
        if(N == 2)
        {
          PyObject * const arg1 = PyTuple_GET_ITEM(_args, 1);
          if(PyInt_Check(arg1)) tolerance = PyInt_AS_LONG(arg1);
          else if(PyFloat_Check(arg1)) tolerance = PyFloat_AS_DOUBLE(arg1);
          else
          {
            LADA_PYERROR(TypeError, "Second argument to cell_invariants should a real number.");
            return NULL;
          }
        }
        return space_group(structure, tolerance);
      }

      PyObject* equivalent_wrapper(PyObject *_module, PyObject *_args, PyObject *_kwargs)
      {
        PyObject *scale = Py_True; 
        PyObject *cartesian = Py_True;
        types::t_real tolerance = types::tolerance;
        PyObject *a = NULL;
        PyObject *b = NULL;
        static char *kwlist[] = { const_cast<char*>("a"),
                                  const_cast<char*>("b"), 
                                  const_cast<char*>("scale"), 
                                  const_cast<char*>("cartesian"), 
                                  const_cast<char*>("tolerance"), NULL};
        if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "OO|OOd:equivalent", kwlist,
                                            &a, &b, &scale, &cartesian, &tolerance) )
          return NULL;
        if(not PyStructure_Check(a))
        {
          LADA_PYERROR(TypeError, "equivalent: First argument should be a structure.");
          return NULL;
        }
        if(not PyStructure_Check(b))
        {
          LADA_PYERROR(TypeError, "equivalent: second argument should be a structure.");
          return NULL;
        }
        if(scale and not PyBool_Check(scale))
        {
          LADA_PYERROR(TypeError, "equivalent: scale should be True or False.");
          return NULL;
        }
        if(cartesian and not PyBool_Check(cartesian))
        {
          LADA_PYERROR(TypeError, "equivalent: cartesian should be True or False.");
          return NULL;
        }
        try
        {
          if( equivalent( Structure::acquire(a), Structure::acquire(b), 
                          scale == Py_True ? true: false,
                          cartesian == Py_True ? true: false,
                          tolerance ) ) Py_RETURN_TRUE;
          Py_RETURN_FALSE;
        }
        // catch exceptions.
        catch(error::pyerror &e)
        {
          // Returns if python error already set.
          if(PyErr_Occurred() != NULL) return NULL;
          // Otherwise, throw our own
          std::string const * message = boost::get_error_info<error::string>(e);
          PyObject ** exception = boost::get_error_info<error::pyexcept>(e);
          if(exception == NULL)
          {
            LADA_PYERROR(InternalError, message ? message->c_str(): 
                             ( "equivalent -- caught following exception: " +
                               std::string(e.what()) ).c_str() );
          }
          else if(message == NULL)
            PyErr_SetString( *exception, 
                             ( "equivalent -- caught following exception: " +
                               std::string(e.what()) ).c_str() );
          else PyErr_SetString(*exception, message->c_str());
        }
        catch(std::exception &e)
        {
          std::string const message =  "equivalent -- caught following exception: " 
                                       + std::string(e.what());
          LADA_PYERROR(InternalError, message.c_str()); 
        }
        catch(...)
        { 
          LADA_PYERROR(InternalError, "Caught unknown c++ exception in equivalent.");
        }
        return NULL;
      }
 
      PyObject* transform_wrapper(PyObject *_module, PyObject *_args, PyObject *_kwargs)
      {
        PyObject *transform_ = NULL; 
        PyObject *structure_ = NULL;
        static char *kwlist[] = { const_cast<char*>("structure"),
                                  const_cast<char*>("transform"), NULL};
        if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "OO:transform", kwlist,
                                            &structure_, &transform_) )
          return NULL;
        if(not PyStructure_Check(structure_))
        {
          LADA_PYERROR(TypeError, "structure: First argument should be a structure.");
          return NULL;
        }
        Eigen::Matrix<types::t_real, 4, 3> transform;
        if(not python::convert_to_matrix(transform_, transform)) return NULL;
        Structure structure = Structure::acquire(structure_).copy();
        structure.transform(transform);
        return structure.release();
     }
      
       
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"zero_centered",  zero_centered_wrapper, METH_VARARGS,
         "Returns Folds vector back to origin.\n\n"
         "This may not be the vector with the smallest possible norm "
           "if the cell is very skewed.\n\n"
         ":Parameters:\n"
         "  a : sequence of three numbers\n"
         "    The vector to fold back into the cell.\n"
         "  cell : sequence of 3x3 numbers\n"
         "    The cell defining the periodicity.\n"
         "  invcell : sequence of 3x3 numbers\n"
         "    Optional. The *inverse* of the cell defining the periodicity. "
              "It is computed if not given on input. \n" },
        {"into_voronoi",  into_voronoi_wrapper, METH_VARARGS,
         "Returns Folds vector into first Brillouin zone of the input cell.\n\n"
         "This returns the periodic image with the smallest possible norm.\n\n"
         ":Parameters:\n"
         "  a : sequence of three numbers\n"
         "    The vector to fold back into the cell.\n"
         "  cell : sequence of 3x3 numbers\n"
         "    The cell defining the periodicity.\n"
         "  invcell : sequence of 3x3 numbers\n"
         "    Optional. The *inverse* of the cell defining the periodicity. "
              "It is computed if not given on input. \n" },
        {"into_cell",  into_cell_wrapper, METH_VARARGS,
         "Returns Folds vector into periodic the input cell.\n\n"
         ":Parameters:\n"
         "  a : sequence of three numbers\n"
         "    The vector to fold back into the cell.\n"
         "  cell : sequence of 3x3 numbers\n"
         "    The cell defining the periodicity.\n"
         "  invcell : sequence of 3x3 numbers\n"
         "    Optional. The *inverse* of the cell defining the periodicity. "
              "It is computed if not given on input. \n" },
        {"are_periodic_images",  are_periodic_images_wrapper, METH_VARARGS,
         "Returns True if first two arguments are periodic images according to third.\n\n"
         ":Parameters:\n"
         "  a : sequence of three numbers\n"
         "    The first vector.\n"
         "  b : sequence of three numbers\n"
         "    The second vector.\n"
         "  invcell : sequence of 3x3 numbers\n"
         "    The *inverse* of the cell defining the periodicity.\n"
         "  tolerance : float\n"
         "    An *optional* floating point defining the tolerance. Defaults to 1e-8." }, 
        {"supercell",  (PyCFunction)supercell_wrapper, METH_VARARGS | METH_KEYWORDS,
         "Creates a supercell of an input lattice.\n\n"
         ":Parameters:\n"
         "  lattice: `crystal.Structure`\n"  
         "    Lattice from which to create the supercell. Cannot be empty. Must be deepcopiable. \n"
         "  cell : 3x3 sequence\n"
         "    Cell in cartesian coordinates of the supercell.\n\n"
         ":returns: A `crystal.Structure` representing the supercell. "
           "It contains an attribute ``lattice`` which reference the input lattice. "
           "If ``lattice`` contains an attribute ``name``, then the result's is set to "
           "\"supercell of ...\". " 
           "All other attributes of the lattice are deep-copied to the supercell. "
           "The atoms within the result contain an attribute ``index`` which is an index to "
           "the equivalent site within ``lattice``. Otherwise, the atoms are deep copies of "
           "the lattice sites. " },
        {"primitive",  primitive_wrapper, METH_VARARGS,
         "Returns the primitive unit-cell lattice of an input lattice.\n\n"
         ":Parameters:\n"
         "  lattice: `crystal.Structure`\n"  
         "    Lattice for which to get primitive unit-cell lattice. Cannot be empty. Must be deepcopiable. \n"
         "  tolerance : float\n"
         "    Optional. Tolerance when comparing positions. Defaults to 1e-8.\n\n"
         ":returns: A new `crystal.Structure` representing the primitive lattice. " },
        {"is_primitive",  is_primitive_wrapper, METH_VARARGS,
         "Returns True if the lattice is primitive.\n\n"
         ":Parameters:\n"
         "  lattice: `crystal.Structure`\n"  
         "    Lattice for which to get primitive unit-cell lattice. Cannot be empty. Must be deepcopiable. \n"
         "  tolerance : float\n"
         "    Optional. Tolerance when comparing positions. Defaults to 1e-8.\n\n"},
        {"space_group",  space_group_wrapper, METH_VARARGS,
         "Finds and stores point group operations.\n\n"
         " Raises a ValueError if the input  structure is not primitive. "
           "Implementation taken from Enum code, PRB 77, 224115 (2008). \n\n"
         ":Parameters:\n"
         "  structure : `lada.crystal.Structure`\n"
         "    The structure for which to find the point group.\n"
         "  tolerance : float\n"
         "    acceptable tolerance when determining symmetries. Defaults to 1e-8.\n"
         ":returns: python list of affine symmetry operations for the given structure. "
           "Each element is a 4x3 numpy array, with the first 3 rows "
           "forming the rotation, and the last row is the translation. "
           "The affine transform is applied as rotation * vector + translation. "},
        {"cell_invariants",  cell_invariants_wrapper, METH_VARARGS,
         "Finds and stores point group operations.\n\n"
         "Rotations are determined from G-vector triplets with the same "
           "norm as the unit-cell vectors. "
           "Implementation taken from Enum code, PRB 77, 224115 (2008). \n\n"
         ":Parameters:\n"
         "  structure : `lada.crystal.Structure`\n"
         "    The structure for which to find the space group. Can also be a 3x3 matrix.\n"
         "  tolerance : float\n"
         "    acceptable tolerance when determining symmetries. Defaults to 1e-8.\n"
         ":returns: python list of affine symmetry operations for the given structure. "
           "Each element is a 4x3 numpy array, with the first 3 rows "
           "forming the rotation, and the last row is the translation. "
           "The affine transform is applied as rotation * vector + translation. "
           "`cell_invariants` always returns rotations (translation is zero). " }, 
        {"equivalent", (PyCFunction)equivalent_wrapper, METH_VARARGS | METH_KEYWORDS, 
          "Returns true if two structures are equivalent. \n\n"
          "Two structures are equivalent in a crystallographic sense, "
          "e.g. without reference to cartesian coordinates or possible "
          "motif rotations which leave the lattice itself invariant. A "
          "supercell is *not* equivalent to its lattice, unless it is a "
          "trivial supercell. An option is provided to allow "
          "comparison within a single reference frame, i.e. per mathematical "
          "definition of lattice.\n\n"
          ":Parameters:\n"
          "  a : `lada.crystal.Structure`\n    The first structure.\n"
          "  b : `lada.crystal.Structure`\n    The second structure.\n"
          "  scale : boolean\n"
          "    whether to take the scale into account. Defaults to true.\n"
          "  cartesian : boolean\n"
          "    whether to take into account differences in cartesian "
              "coordinates. Defaults to true. If False, then comparison is "
              "according to mathematical definition of a lattice. If True, "
              "comparison is according to crystallographic comparison.\n"
          "  tolerance : float\n"
          "    Tolerance when comparing distances. Defaults to 1e-8.  It is in "
              "the same units as the structures scales, if that is taken into "
              "account, otherwise, it is in the same units as ``a.scale``." },
        {"transform", (PyCFunction)transform_wrapper, METH_VARARGS | METH_KEYWORDS, 
          "Returns a copy of the structure transformed according to affine operation. \n\n"
          ":Parameters:\n"
          "  structure : `lada.crystal.Structure`\n    The structure to transform.\n"
          "  transform : 4x3 matrix\n    The rotation (applied first) is contained "
            "in the first three row and the translation (applied second) in the last row." },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    }; // end of static method table.
  }
}

