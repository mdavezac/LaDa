namespace LaDa
{
  namespace crystal
  {
    //! Wrapper around the supercell method.
    PyObject* supercell_wrapper(PyObject *_lattice, PyObject *_cell)
    {
      // check/convert input parameters.
      if(not PyStructure_Check(_lattice))
      {
        LADA_PYERROR(TypeError, "Input is not a crystal.Structure object."); 
        return NULL;
      }
      Structure lattice((StructureData*)_lattice);
      math::rMatrix3d cell;
      if(not math::convert_to_cell(_cell, cell)) return NULL;
      // create supercell.
      try { return supercell((StructureData*)_lattice, cell).release(); } 
      // catch exceptions.
      catch(error::ValueError &e)
      {
        std::string const * error = boost::get_error_info<error::string>(e);
        LADA_PYERROR(ValueError, error ? error->c_str(): "");
      }
      catch(error::internal &e)
      {
        std::string const * error = boost::get_error_info<error::string>(e);
        LADA_PYERROR(ValueError, error ? error->c_str(): "");
      }
      catch(std::exception &e) { LADA_PYERROR(internal, e.what()); }
      catch(...) { LADA_PYERROR(internal, "Caught unknown c++ exception."); }
      return NULL;
    }
       
    //! Methods table for crystal module.
    static PyMethodDef methods_table[] = {
        {"supercell",  supercell_wrapper, METH_VARARGS,
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
           "All other attributes of the lattice are copied to the supercell. "
           "The atoms within the result contain an attribute ``index`` which is an index to "
           "the equivalent site within ``lattice``. Otherwise, the atoms are deep copies of "
           "the lattice sites. " },
        {NULL, NULL, 0, NULL}        /* Sentinel */
    }; // end of static method table.
  }
}

