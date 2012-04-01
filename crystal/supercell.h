#ifndef LADA_CRYSTAL_SUPERCELL_H
#define LADA_CRYSTAL_SUPERCELL_H
#include "LaDaConfig.h"

#include <cmath>

#include <Eigen/LU> 

#include <math/misc.h>
#include <math/fuzzy.h>

#include "smith/smith.h"
#include "utilities.h"


namespace LaDa
{
  namespace crystal 
  {

    template<class T_DERIVED>
      Structure supercell(Structure const &_lattice, Eigen::DenseBase<T_DERIVED> const &_supercell)
      {
        if(_lattice.size() == 0) 
        {
          LADA_PYTHROW(ValueError, "Lattice is empty.");
          return Structure();
        }
        namespace bt = boost::tuples;
        Structure result = _lattice.copy(); 
        if(not result) LADA_PYTHROW(ValueError, "Could not deepcopy the lattice.");
        result.clear();
        result->cell = _supercell;
        if(_lattice.hasattr("name") and PyString_Check(_lattice.pyattr("name").borrowed()) )
        {
          char *const attr = PyString_AS_STRING(_lattice.pyattr("name").borrowed());
          if(attr != "" and not result.pyattr_convert("name", "supercell of " + std::string(attr)) ) 
            PyErr_Clear();
        }
        SmithTransform transform( _lattice.cell(), result.cell());
        if(not transform) return Structure();;
      
        const math::rMatrix3d factor(transform.transform().inverse());
        math::rMatrix3d inv_cell( result.cell().inverse() ); 
        result.reserve(transform.size()*_lattice.size());
        Structure::const_iterator const i_site_begin = _lattice.begin();
        Structure::const_iterator const i_site_end = _lattice.end();
        
        for( math::iVector3d::Scalar i(0); i < transform.quotient()(0); ++i )
          for( math::iVector3d::Scalar j(0); j < transform.quotient()(1); ++j )
            for( math::iVector3d::Scalar k(0); k < transform.quotient()(2); ++k )
            {
              // in cartesian.
              const math::rVector3d vec( factor * math::rVector3d(i,j,k) );
            
              // adds all lattice sites.
              long l(0);
              for(Structure::const_iterator i_site(i_site_begin); i_site != i_site_end; ++i_site, ++l)
              {
                Atom atom = i_site->copy();
                if(not atom) LADA_PYTHROW(ValueError, "Could not deepcopy atom.");
                atom->pos = into_cell(vec+i_site->pos(), result->cell, inv_cell);
                python::Object site = PyInt_FromLong(l);
                if(not atom) LADA_PYTHROW(internal, "Could not create python integer.");
                if(not atom.pyattr("site", site) ) 
                  LADA_PYTHROW(internal, "Could not set site index attribute.");
                result.push_back(atom);
              }
            }
      
        return result;
      }
  } // namespace crystal
} // namespace LaDa
#endif
