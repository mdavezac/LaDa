#ifndef LADA_CRYSTAL_MAPSITES_H
#define LADA_CRYSTAL_MAPSITES_H

#include "LaDaConfig.h"

#include <math/misc.h>
#include <math/gruber.h>

#include "utilities.h"
#include "structure.h"
#include "exceptions.h"
#include "compare_sites.h"


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
        if(_mapper.size() == 0) 
          BOOST_THROW_EXCEPTION(error::input() << error::string("Empty mapper structure."));
        if(_mappee.size() == 0) 
          BOOST_THROW_EXCEPTION(error::input() << error::string("Empty mappee structure."));

        math::rMatrix3d const cell = math::gruber(_mapper.cell());
        math::rMatrix3d const invcell = cell.inverse();
        
        // check that mappee_ is a supercell of mapper_.
        types::t_real const ratio = _mappee->scale / _mapper->scale;
        types::t_real tolerance = _tolerance / _mapper->scale;
        math::rMatrix3d const intcell_ = invcell * _mappee.cell() * ratio;
        if(not math::is_integer(intcell_, _tolerance))
          BOOST_THROW_EXCEPTION(error::not_a_supercell());

        // Copy mapper sites to a vector, making sure positiosn are in cell.
        std::vector<math::rVector3d> sites; 
        typename TemplateStructure<T_TYPE>::const_iterator i_mapper_site = _mapper.begin();
        typename TemplateStructure<T_TYPE>::const_iterator const i_mapper_site_end = _mapper.end();
        for(; i_mapper_site != i_mapper_site_end; ++i_mapper_site)
          sites.push_back(into_cell(i_mapper_site->pos(), cell, invcell));

        // loop over atoms in mappee and assign sites.
        bool allmapped = true;
        typename TemplateStructure<T_TYPE>::iterator i_atom = _mappee.begin();
        typename TemplateStructure<T_TYPE>::iterator const i_atom_end = _mappee.end();
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
              = into_voronoi(ratio*i_atom->pos()-(*i_site), cell, invcell).squaredNorm();
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
            BOOST_THROW_EXCEPTION(error::two_sites_at_same_position());
          if(fneigh_dist > tolerance) { i_atom->site() = -1; allmapped = false; }
          else if(_withocc and (not compare_occupations(_mapper[fneigh_index]->type)(i_atom->type())))
            { i_atom->site() = -1; allmapped = false; }
          else i_atom->site() = fneigh_index;
        }
        return allmapped;
      }
  } // namespace crystal
} // namespace LaDa

#endif
