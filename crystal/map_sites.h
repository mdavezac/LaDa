#ifndef LADA_CRYSTAL_MAPSITES_H
#define LADA_CRYSTAL_MAPSITES_H

#include "LaDaConfig.h"

#include <math/misc.h>
#include <math/gruber.h>

#include "structure.h"
#include "neighbors.h"
#include "exceptions.h"
#include "supercell.h"


namespace LaDa 
{
  namespace crystal
  {
    //! \brief Map atomic sites from mapper onto mappee.
    //! \param[in] _mapper : a lattice against which to map atomic sites.
    //! \param[inout] _mappee : a supercell for which sites will be mapped.
    //! \param[in] _withocc : whether to take occupation into count, or only position.
    //! \param[in] _tolerance : Tolerance criteria for distances, in units of _mappee.scale().
    //! \details Where possible, the site indices of the mappee structure
    //!          corresponds to the equivalent sites in the mapper structure.
    //! \return True if mapping is successful, False if all sites could not be mapped. 
    //!         Since in the case of defects, incomplete mappings may be what is wanted, 
    template<class T_TYPE>
      bool map_sites( TemplateStructure<T_TYPE> const &_mapper, TemplateStructure<T_TYPE> &_mappee,
                      bool _withocc = true, types::t_real _tolerance = types::tolerance )
      {
        math::rMatrix3d const cell = math::gruber(_mapper.cell());
        math::rMatrix3d const intcell_ = cell.inverse() * _mappee.cell();
        if(not math::is_integer(intcell_, _tolerance))
          BOOST_THROW_EXCEPTION(error::not_a_supercell());
        math::rMatrix3d intcell;
        for(size_t i(0); i < 3; ++i)
          for(size_t j(0); j < 3; ++j)
            intcell(i, j) = std::floor(intcell_(i,j) + 0.2);
        TemplateStructure<T_TYPE> lattice = supercell(_mapper, cell * intcell);
        lattice.cell() *= lattice.scale();
        lattice.scale() = 1e0;
        math::rMatrix3d const transform = lattice.cell() * _mappee.cell().inverse() * _mappee.scale();

        types::t_real tolerance = _tolerance * _mappee.scale();
        bool allmapped = true;
        typename TemplateStructure<T_TYPE>::iterator i_atom = _mappee.begin();
        typename TemplateStructure<T_TYPE>::iterator const i_atom_end = _mappee.end();
        for(; i_atom != i_atom_end; ++i_atom)
        {
          Neighbors neighbors(2, transform * i_atom->pos, _tolerance, true);
          Neighbors::const_iterator i_first = neighbors.begin(lattice);
          Neighbor const & first = *i_first;
          Neighbor const & second = *(++i_first);
          std::cout << first.distance << " " << second.distance << "\n";
          if( math::eq(first.distance, second.distance, tolerance) )
            BOOST_THROW_EXCEPTION(error::two_sites_at_same_position());
          if(first.distance > tolerance)
            { i_atom->site = -1; allmapped = false; }
          else if(_withocc and (not compare_occupations(lattice[first.index].type)(i_atom->type)))
            { i_atom->site = -1; allmapped = false; }
          else i_atom->site = lattice[first.index].site;
        }
        return allmapped;
      }
  } // namespace crystal
} // namespace LaDa

#endif
