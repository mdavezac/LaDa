//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iterator>

#include "variable_major.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {

      CoordRange::rank_range
        ::rank_range( SumOfSeparables::const_iterator const &_irank,
                      SumOfSeparables::const_iterator const &_irank_end,
                      size_t _i)
            : i_rank_(_irank), i_rank_end_(_irank_end), index_(_i) 
      {
        try
        {
          while( i_rank_->get<0>().size() <= index_  )
          {
            ++i_rank_;
            if( i_rank_ == i_rank_end_ ) break;
          } 
        }
        catch(...)
        {
          i_rank_ == i_rank_end_;
        }
      }

      CoordRange::rank_range::var_iterator_ 
        CoordRange::rank_range::get_i_var_() const 
        {
          LADA_ASSERT( i_rank_ != i_rank_end_, "Range not iteratable.\n" )
          var_iterator_ result( i_rank_->get<0>().begin() );
          for(size_t i(0); i < index_; ++i, ++result);
          return result;
        };

      void CoordRange::rank_range::operator++() 
      {
        if( i_rank_ != i_rank_end_ ) return;

        do
        {
          ++i_rank_;
          if( i_rank_ == i_rank_end_ ) break;

        } while( i_rank_->get<0>().size() <= index_  );

      }

    } // namespace collapse.
  } // namespace atomic_potential.
} // namespace LaDa
