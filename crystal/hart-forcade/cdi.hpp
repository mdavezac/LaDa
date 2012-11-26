//! Function to initialize a string atom.
int hftransform_init(PyHFTObject* _self, PyObject* _args, PyObject *_kwargs);


//! \brief Initializes a new hftransform from input lattice unit-cell and supercell.
//! \details Performs initialization from c++ arguments.
bool _init_hft( PyHFTObject* _self, 
                           math::rMatrix3d const &_lattice,
                           math::rMatrix3d const &_supercell )
{
  if(std::abs(_lattice.determinant()) < 1e-8)
  {
    LADA_PYERROR(ValueError, "Unit-cell is singular.");
    return false;
  }
  if(std::abs(_supercell.determinant()) < 1e-8)
  {
    LADA_PYERROR(ValueError, "Supercell is singular.");
    return false;
  }
  math::iMatrix3d left, right, hf;
  const math::rMatrix3d inv_lat( !_lattice );
  const math::rMatrix3d inv_lat_cell( inv_lat * _supercell );
  math::iMatrix3d int_cell;
  for( size_t i(0); i < 3; ++i )
    for( size_t j(0); j < 3; ++j )
    {
      int_cell(i,j) = types::t_int( rint( inv_lat_cell(i,j) ) ); 
      if( math::neq(types::t_real(int_cell(i,j)), inv_lat_cell(i,j), 1e-2) )
      {
        LADA_PYERROR(ValueError, "Second argument is not a supercell of first.");
        return false;
      }
    }
  try
  { 
    math::smith_normal_form( hf, left, int_cell, right );
    _self->transform = left.cast<types::t_real>() * (!_lattice);
    _self->quotient = hf.diagonal();
    return true;
  }
  catch(error::internal &_e)
  {
    LADA_PYERROR(internal, "HFTransform: Could not create smith normal form.");
  }
  catch(std::exception &_e)
  {
    LADA_PYERROR_FORMAT(internal, "HFTransform: Caught c++ exception %s.", _e.what());
  }
  catch(...)
  {
    LADA_PYERROR(internal, "HFTransform: Caught unknown c++ exception.");
  }
  return false;
}

// Function to initialize an atom.
int hftransform_init(PyHFTObject* _self, PyObject* _args, PyObject *_kwargs)
{
  PyObject *lattice = NULL;
  PyObject *supercell = NULL;
  static char *kwlist[] = { const_cast<char*>("unitcell"), const_cast<char*>("supercell"), NULL };
  if(not PyArg_ParseTupleAndKeywords( _args, _kwargs, "OO:HFTransfrom", kwlist,
                                      &lattice, &supercell ) )
    return -1;
  math::rMatrix3d cell, bigcell;
  if(check_structure(lattice)) cell = ((PyStructureObject*)lattice)->cell;
  else if(not python::numpy::convert_to_matrix(lattice, cell)) { std::cout << "Err 0\n"; return -1; }
  if(check_structure(supercell)) bigcell = ((PyStructureObject*)supercell)->cell;
  else if(not python::numpy::convert_to_matrix(supercell, bigcell)){ std::cout << "Err 1\n"; return -1; }
  
  return _init_hft(_self, cell, bigcell) ? 0: -1;
}


