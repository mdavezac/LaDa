#ifndef LADA_CRYSTAL_SPACEGROUP_H
#define LADA_CRYSTAL_SPACEGROUP_H

#include "LaDaConfig.h"

#include <vector>
#include <iostream>

#include <boost/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>
#include <Eigen/StdVector>

#include <opt/types.h>
#include <math/fuzzy.h>
#include <math/misc.h>
#include <math/gruber.h>

#include "structure.h"
#include "exceptions.h"
#include "primitive.h"

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

    //! \brief Finds and stores space group operations.
    //! \param[in] _structure The structure for which to find the space group.
    //! \param[in] _tol acceptable tolerance when determining symmetries.
    //!             -1 implies that types::tolerance is used.
    //! \retval spacegroup Shared pointer to a t_SpaceGroup vector containing
    //!         the symmetry operations for the given structure.
    //! \warning Works for primitive lattices only.
    //! \see Taken from Enum code, PRB 77, 224115 (2008).
    template<class T_TYPE> boost::shared_ptr<t_SpaceGroup>
      space_group(TemplateStructure<T_TYPE> const &_lattice, types::t_real _tolerance = -1e0)
      {
        if(_tolerance <= 0e0) _tolerance = types::tolerance;
        // Checks that lattice has atoms.
        if(_lattice.size() == 0) BOOST_THROW_EXCEPTION(error::empty_structure());
        // Checks that lattice is primitive.
        if(not _lattice.is_clone(primitive(_lattice))) BOOST_THROW_EXCEPTION(error::not_primitive());
        
 
        // Finds minimum translation.
        TemplateStructure<T_TYPE> atoms = _lattice.copy();
        math::rVector3d translation(atoms.front().pos);
        math::rMatrix3d const cell(math::gruber(_lattice.cell()));
        math::rMatrix3d const invcell(!cell);
        // Creates a list of atoms centered in the cell.
        foreach(typename TemplateStructure<T_TYPE>::reference site, atoms)
          site.pos = into_cell(site.pos-translation, cell, invcell);
 
        // gets point group.
        boost::shared_ptr<t_SpaceGroup> pg = cell_invariants(_lattice.cell());
        boost::shared_ptr<t_SpaceGroup> result(new t_SpaceGroup);
        result->reserve(pg->size());
             
        // lists atoms of same type as atoms.front()
        std::vector<math::rVector3d> translations;
        CompareOccupations<T_TYPE> const compsites(atoms.front().type);
        foreach(typename TemplateStructure<T_TYPE>::const_reference site, atoms)
          if(compsites(site.type)) translations.push_back(site.pos);
        
 
        // applies point group symmetries and finds out if they are part of the space-group.
        foreach(t_SpaceGroup::reference op, *pg)
        {
          // loop over possible translations.
          std::vector<math::rVector3d> :: const_iterator i_trial = translations.begin();
          std::vector<math::rVector3d> :: const_iterator const i_trial_end = translations.end();
          typedef typename TemplateStructure<T_TYPE>::const_iterator const_iterator;
          const_iterator const i_atoms_begin = atoms.begin();
          const_iterator const i_atoms_end = atoms.end();
          for(; i_trial != i_trial_end; ++i_trial)
          {
            // possible translation.
            op.translation() = *i_trial;
            // Checks that this is a mapping of the lattice upon itself.
            const_iterator i_unmapped = i_atoms_begin;
            for(; i_unmapped != i_atoms_end; ++i_unmapped)
            {
              CompareSites<T_TYPE> const transformed( into_cell(op*i_unmapped->pos, cell, invcell), 
                                                      i_unmapped->type, _tolerance );
              const_iterator i_mapping = i_atoms_begin;
              for(; i_mapping != i_atoms_end; ++i_mapping)
                if(transformed(*i_mapping)) break;
              // found unmapped site if condition is true.
              if(i_mapping == i_atoms_end) break;
            } // loop over all atoms.
 
            // all sites in the lattice were mapped if condition is true.
            if(i_unmapped == i_atoms_end) break; 
          } // loop over trial translations.
 
          // Found transformation which maps lattice upon itself if condition is true.
          if(i_trial != i_trial_end) result->push_back(op);
        } // loop over point group.
 
        return result;
      } 
  }
}

#endif
