#ifndef LADA_CRYSTAL_SYMMETRY_OPERATOR_H_
#define LADA_CRYSTAL_SYMMETRY_OPERATOR_H_

#include "LaDaConfig.h"

#include <vector>
#include <iostream>

#include <boost/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>
#include <Eigen/StdVector>

#include <opt/debug.h>
#include <opt/types.h>
#include <math/fuzzy.h>
#include <math/misc.h>
#include <math/serialize.h>


namespace LaDa
{
  namespace crystal 
  {
    //! \typedef Vector of symmetry operations making up a space-group.
    typedef std::vector<math::Affine3d, Eigen::aligned_allocator<math::Affine3d> > t_SpaceGroup;

    //! \brief Returns point symmetries of a cell (except identity).
    //! \details Rotations are determined from G-vector triplets with the same
    //!          norm as the unit-cell vectors.
    //! \see Taken from Enum code, PRB 77, 224115 (2008).
    boost::shared_ptr<t_SpaceGroup>
      cell_invariants(math::rMatrix3d const &_cell, types::t_real _tolerance = -1e0);

//   //! \brief Finds and stores space group operations.
//   //! \param[in] _structure The structure for which to find the space group.
//   //! \param[in] _tol acceptable tolerance when determining symmetries.
//   //!             -1 implies that types::tolerance is used.
//   //! \retval spacegroup Shared pointer to a t_SpaceGroup vector containing
//   //!         the symmetry operations for the given structure.
//   //! \warning Works for primitive lattices only.
//   //! \see Taken from Enum code, PRB 77, 224115 (2008).
//   template<class T_TYPE> boost::shared_ptr<t_SpaceGroup>
//     space_group(TemplateStructure<T_TYPE> const &_lattice, types::t_real _tolerance = -1e0);
//     {
//       if(_tolerance <= 0e0) _tolerance = types::tolerance;
//       // Checks that lattice has atoms.
//       if(_lattice.size() == 0) BOOST_THROW_EXCEPTION(error::empty_structure());
//       { // Checks that lattice is primitive.
//         Lattice lat(_lattice);
//         LADA_DOASSERT( lat.make_primitive(), 
//                        "Lattice is not primitive.\nCannot compute symmetries.\n" );
//       }
//
//       // Finds minimum translation.
//       Lattice::t_Sites atoms(_lattice.atoms);
//       math::rVector3d translation(atoms.front().pos);
//       math::rMatrix3d const invcell(!_lattice.cell);
//       // Creates a list of atoms centered in the cell.
//       foreach( Lattice::t_Site &site, atoms )
//         site.pos = into_cell(site.pos-translation, _lattice.cell, invcell);
//
//       // gets point group.
//       boost::shared_ptr< std::vector<SymmetryOperator> >
//         pg = get_point_group_symmetries(_lattice.cell);
//       boost::shared_ptr< std::vector<SymmetryOperator> > 
//         result( new std::vector<SymmetryOperator> );
//       result->reserve(pg->size());
//            
//       // lists atoms of same type as atoms.front()
//       std::vector<math::rVector3d> translations;
//       CompareSites compsites(atoms.front(), _tolerance);
//       foreach( Lattice::t_Site const &site, atoms )
//         if( compsites(site.type) ) translations.push_back(site.pos);
//       
//
//       // applies point group symmetries and finds out if they are part of the space-group.
//       foreach(SymmetryOperator &op, *pg)
//       {
//         // loop over possible translations.
//         std::vector<math::rVector3d> :: const_iterator i_trial = translations.begin();
//         std::vector<math::rVector3d> :: const_iterator const i_trial_end = translations.end();
//         for(; i_trial != i_trial_end; ++i_trial)
//         {
//           // possible translation.
//           op.trans = *i_trial;
//           Lattice::t_Sites::const_iterator i_site = atoms.begin();
//           Lattice::t_Sites::const_iterator const i_site_end = atoms.end();
//           for(; i_site != i_site_end; ++i_site)
//           {
//             CompareSites transformed(*i_site, _tolerance);
//             transformed.pos = into_cell(op(i_site->pos), _lattice.cell, invcell);
//             Lattice::t_Sites::const_iterator const
//               i_found( std::find_if(atoms.begin(), atoms.end(), transformed) );
//             if(i_found == atoms.end()) break;
//             if( not transformed(i_found->type) ) break;
//           } // loop over all atoms.
//
//           if(i_site == i_site_end) break; // found a mapping
//         } // loop over trial translations.
//
//         if(i_trial != i_trial_end)
//           result->push_back(SymmetryOperator(op.op, op.trans-op.op*translation+translation) );
//       } // loop over point group.
//
//       return result;
//     } 
  }
}

#endif
