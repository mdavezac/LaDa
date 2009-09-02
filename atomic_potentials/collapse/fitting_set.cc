//
//  Version: $Id: collapse.cc 1270 2009-08-17 03:11:58Z davezac $
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../representation.h"
#include "fitting_set.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {
      FittingSet::str_iterator FittingSet::begin(size_t _i) const
      {
        LADA_ASSERT( _i < coordinates_.size(), "Index out of range.\n" );

        str_iterator result;
        result.i_ = _i;
        result.i_coordinates_ = coordinates_[_i].begin();
        result.i_energy_ = energies_.begin();
        result.i_weight_ = weights_.begin();
        return result;
      }
      FittingSet::str_iterator FittingSet::end(size_t _i) const
      {
        LADA_ASSERT( _i < coordinates_.size(), "Index out of range.\n" );

        str_iterator result;
        result.i_ = _i;
        result.i_coordinates_ = coordinates_[_i].end();
        result.i_energy_ = energies_.end();
        result.i_weight_ = weights_.end();
        return result;
      }

      FittingSet::str_iterator::rep_iterator
        FittingSet::str_iterator::begin() const
        {
          rep_iterator result;
          result.i_coordinates_ = i_coordinates_->begin();
          result.i_weight_ = i_weight_->second.begin();
          return result;
        }

      FittingSet::str_iterator::rep_iterator
        FittingSet::str_iterator::end() const
        {
          rep_iterator result;
          result.i_coordinates_ = i_coordinates_->end();
          result.i_weight_ = i_weight_->second.end();
          return result;
        }

      void FittingSet::add( Representation const &_representation,
                            numeric_type _energy, numeric_type _weight )
      {
        LADA_ASSERT( _representation.size(), "Empty representation.\n")
        LADA_ASSERT( _representation.nb_coords(), "Empty representation.\n")
        // adds types.
        if( coordinates_.size() == 0 ) coordinates_.resize( _representation.nb_coords() );
        t_Coordinates::iterator i_coord = coordinates_.begin();
        t_Coordinates::iterator const i_coord_end = coordinates_.end();
        for(; i_coord != i_coord_end; ++i_coord)
        { 
          i_coord->resize( i_coord->size()+1 );
          t_Coordinates::value_type::value_type& coord_container = i_coord->back();
          Representation::const_iterator i_rep = _representation.begin();
          Representation::const_iterator const i_rep_end = _representation.end();
          for(; i_rep != i_rep_end; ++i_rep)
          {
            // creates coordinate vector.
            i_coord->back().resize(i_coord->back().size()+1);
            t_Coordinates::value_type::value_type
                                     ::value_type &type_container = i_coord->back().back();
            foreach(VariableSet::t_Variables::value_type const& value, i_rep->variables)
              type_container.push_back( value.second );
          }
        }
        // adds energy.
        energies_.push_back(_energy);

        // adds weights.
        {
          Representation::const_iterator i_rep = _representation.begin();
          Representation::const_iterator const i_rep_end = _representation.end();
          t_Weights::value_type::second_type weights;
          for(; i_rep != i_rep_end; ++i_rep) weights.push_back(i_rep->weight);
          weights_.push_back( t_Weights::value_type(_weight, weights) );
        }
      }
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
