#include "LaDaConfig.h"

#include <math/misc.h>
#include <math/gruber.h>
#include "utilities.h"
#include "map_sites.h"

namespace LaDa 
{
  namespace crystal
  {

    bool map_sites( Structure const &_mapper, Structure &_mappee,
                    python::Object _withocc, types::t_real _tolerance )
    {
      if(_mapper.size() == 0) 
      {
        LADA_PYERROR(ValueError, "Empty mapper structure.");
        return false;
      }
      if(_mappee.size() == 0) 
      {
        LADA_PYERROR(ValueError, "Empty mappee structure.");
        return false;
      }

      math::rMatrix3d const cell = math::gruber(_mapper.cell());
      math::rMatrix3d const invcell = cell.inverse();
      bool withocc( PyCallable_Check(_withocc.borrowed()) );
      
      // check that mappee_ is a supercell of mapper_.
      types::t_real const ratio = _mappee->scale / _mapper->scale;
      types::t_real tolerance = _tolerance / _mapper->scale;
      math::rMatrix3d const intcell_ = invcell * _mappee.cell() * ratio;
      if(not math::is_integer(intcell_, _tolerance))
      {
        LADA_PYERROR(ValueError, "Mappee not a supercell of mapper.");
        return false;
      }

      // Copy mapper sites to a vector, making sure positiosn are in cell.
      std::vector<math::rVector3d> sites; 
      Structure::const_iterator i_mapper_site = _mapper.begin();
      Structure::const_iterator const i_mapper_site_end = _mapper.end();
      for(; i_mapper_site != i_mapper_site_end; ++i_mapper_site)
        sites.push_back(into_cell(i_mapper_site->pos(), cell, invcell));

      // loop over atoms in mappee and assign sites.
      bool allmapped = true;
      Structure::iterator i_atom = _mappee.begin();
      Structure::iterator const i_atom_end = _mappee.end();
      std::vector<math::rVector3d>::const_iterator const i_site_end = sites.end();
      for(; i_atom != i_atom_end; ++i_atom)
      {
        // loop over lattice sites, find two first neighbors.
        types::t_int fneigh_index = -1;
        types::t_int sneigh_index = -1;
        types::t_real fneigh_dist = -1;
        types::t_real sneigh_dist = -1;
        std::vector<math::rVector3d>::const_iterator i_site = sites.begin();
        for(size_t i(0); i_site != i_site_end; ++i_site, ++i)
        {
          types::t_real const norm 
            = math::absnorm(into_voronoi(ratio*i_atom->pos()-(*i_site), cell, invcell));
          if(fneigh_dist > norm or fneigh_index == -1) 
          {
            sneigh_dist = fneigh_dist;
            sneigh_index = fneigh_index;
            fneigh_dist = norm;
            fneigh_index = i;
          }
          else if(sneigh_dist > norm or sneigh_index == -1)
          {
            sneigh_dist = norm;
            sneigh_index = i;
          }
        }
        if( math::eq(fneigh_dist, sneigh_dist, tolerance) and sneigh_index != -1)
        {
          LADA_PYERROR(ValueError, "Found two atoms at the same site.");
          return false;
        }
        if(fneigh_dist > tolerance) 
        {
          fneigh_index = -1;
          allmapped = false; 
        }
        else if(withocc) 
        {
          python::Object result
            = PyObject_CallFunctionObjArgs( _withocc.borrowed(),
                                            _mapper[fneigh_index].borrowed(),
                                            i_atom->borrowed(), NULL );
          if(not result) BOOST_THROW_EXCEPTION(error::internal());
          if(PyBool_Check(result.borrowed()))
          {
            if(result.borrowed() == Py_False) fneigh_index = -1;
          }
          else if(PyLong_Check(result.borrowed()) or PyInt_Check(result.borrowed()))
          {
            int i = PyInt_AS_LONG(result.borrowed());
            if(i == -1 and PyErr_Occurred() != NULL) BOOST_THROW_EXCEPTION(error::internal());
            if(i == 0) fneigh_index = -1;
          }
          else
          {
            LADA_PYERROR(ValueError, "Callable is expected to return True or False");
            return false;
          }
        }
        if(fneigh_index == -1)
        {
          i_atom->pyattr("site", Py_None);
          if(PyErr_Occurred() != NULL) BOOST_THROW_EXCEPTION(error::internal());
        }
        else
        { 
          python::Object pyint = PyLong_FromLong(fneigh_index);
          if(not pyint)
          { 
            PyErr_Clear(); 
            LADA_PYERROR(internal, "Could not create python integer."); 
            return false;
          }
          i_atom->pyattr("site", pyint.borrowed());
        }
      }
      return allmapped;
    }
  } // namespace crystal
} // namespace LaDa

