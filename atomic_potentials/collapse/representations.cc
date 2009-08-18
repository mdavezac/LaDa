//
//  Version: $Id: collapse.cc 1270 2009-08-17 03:11:58Z davezac $
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "representation.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {
      Representations::str_iterator Representations::begin(size_t _i) const
      {
        LADA_ASSERT( _i < coordinates_.size(), "Index out of range.\n" );

        str_iterator result;
        result.i_coordinates_ = coordinates_[_i].begin();
        result.i_energy = energies_.begin();
        result.i_weight_ = weights_.begin();
        return result;
      }
      Representations::str_iterator Representations::end(size_t _i) const
      {
        LADA_ASSERT( _i < coordinates_.size(), "Index out of range.\n" );

        str_iterator result;
        result.i_coordinates_ = coordinates_[_i].end();
        result.i_energy = energies_.end();
        result.i_weight_ = weights_.end();
        return result;
      }

      Representations::str_iterator::rep_iterator
        Representations::str_iterator::begin() const
        {
          rep_iterator result;
          result.i_coordinates_ = i_coordinates_->begin();
          result.i_weight_ = weights_->second.begin();
          return result;
        }

      Representations::str_iterator::rep_iterator
        Representations::str_iterator::end() const
        {
          rep_iterator result;
          result.i_coordinates_ = i_coordinates_->end();
          result.i_weight_ = weights_->second.end();
          return result;
        }
    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
