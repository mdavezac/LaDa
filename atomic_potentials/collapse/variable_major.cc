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

      CoordRange::CoordRange(SumOfSeparables const &_sumofseps ) 
         : sumofseps_(_sumofseps), i_(0), max_index_(0) {}
      {
        typedef SumOfSeparables::t_Functions::const_iterator const_iterator;
        const_iterator i_sep( sumofseps_.functions_.begin() );
        const_iterator const i_sep_end( sumofseps_.functions_.end() );
        for(; i_func != i_func_end; ++i_func)
          max_index_ = std::max( max_index_, i_func->size() );
      }  
      
      CoordRange::rank_range
        ::rank_iterator( SumOfSeparables::const_iterator const &_irank,
                         SumOfSeparables::const_iterator const &_irank_end,
                         size_t _i)
            : i_rank_(_irank), i_rank_end_(_irank_end), index_(_i) {}
      {
        try
        {
          while( i_rank->get<0>().size() <= index_  )
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
        CoordRange::rank_range::get_i_var_() const;
        {
          LADA_ASSERT( i_rank_ != i_rank_end, "Range not iteratable.\n" )
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

        } while( i_rank->get<0>().size() <= index_  );

      }

    } // namespace collapse.
  } // namespace atomic_potential.
} // namespace LaDa
