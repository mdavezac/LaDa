#ifndef _LADA_CRYSTAL_SMITH_H_
#define _LADA_CRYSTAL_SMITH_H_

#include "LaDaConfig.h"

#include <boost/tuple/tuple.hpp>

#include <math/smith_normal_form.h>
#include "structure.h"

namespace LaDa 
{
  namespace crystal 
  {
    //! Computes linear smith index of position \a _pos.
    template<class T_TYPE>
      size_t linear_smith_index( t_SmithTransform const &_transformation,
                                 Crystal::Atom_Type<T_TYPE> &_atom ) 
      {
        if(_atom.site < 0) BOOST_THROW_EXCEPTION( error::incorrect_site_index() );
        return math::linear_smith_index
               (
                 boost::tuples::get<1>(_transformation),
                 _atom.site,
                 get_smith_index( _transformation, _atom.pos )
               );
      }


    //! \brief Computes map of smith indices.
    //! \param [in] _lattice: Backbone lattice.
    //! \param [in] _supercell: Supercell of the lattice. The site indices of
    //!                         each atom must correspond to the index of the
    //!                         site in the backbone lattice.
    //! \param [out] _map: Two-dimensional matrix where rows (inner vector) corresponds to
    //!                    atomic sites indexed by their linear smith index,
    //!                    and columns (outer vector) to lattice sites. The
    //!                    value of each element is the corresponding index
    //!                    into the supercell's list of atoms.
    template<class T_TYPE>
      void smith_map( Crystal::TemplateStructure<T_TYPE> const &_lattice,
                      Crystal::TemplateStructure<T_TYPE> const &_supercell,
                      std::vector< std::vector<size_t> > &_map )
      {
        // finds map for atomic positions.
        Crystal::t_SmithTransform const
          transform( get_smith_transform(_lattice.cell(), _supercell.cell()) );
     
        if(_supercell.size() % _lattice.size() != 0)
          BOOST_THROW_EXCEPTION(error::incommensurate_number_of_sites());
        size_t const Nat( _str.atoms.size() );
        size_t const N( Nat / _str.lattice->sites.size() );
        _map.clear();
        _map.resize(_str.lattice->sites.size(), std::vector<size_t>(N, Nat)); 
     
        Crystal::Structure::t_Atoms::const_iterator i_first = _str.atoms.begin();
        Crystal::Structure::t_Atoms::const_iterator const i_end = _str.atoms.end();
        for(size_t i(0); i_first != i_end; ++i_first, ++i)
        {
          if(i_first->site < 0 or i_first->site >= _str.lattice->sites.size())
            BOOST_THROW_EXCEPTION(incorrect_site_index());
          
          math::rVector3d const pos( i_first->pos - _str.lattice->sites[i_first->site].pos );
          size_t const smith( Crystal::get_linear_smith_index(transform, pos) );
#         ifdef LADA_DEBUG
            if(smith >= N) BOOST_THROW_EXCEPTION(error::internal());
            if(_map[i_first->site][smith]!=Nat) BOOST_THROW_EXCEPTION(error::internal());
#         endif
          _map[i_first->site][smith] = i;
        }
#       ifdef LADA_DEBUG
          for( size_t site(0), j(0); site < _str.lattice->sites.size(); ++site )
            for( size_t i(0); i < N; ++i, ++j )
              if(_map[site][i] == Nat) BOOST_THROW_EXCEPTION(error::internal())
#       endif
    }
  } // namespace Crystal

} // namespace LaDa

#endif
