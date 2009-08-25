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

#     ifdef LADA_GETIVAR
#       error LADA_GETIVAR already defined.
#     endif
#     define LADA_GETIVAR( a, b, c)  \
        CoordRange::a ## rank_range::b ## var_iterator \
          CoordRange::rank_range::b ## get_i_var() c \
          { \
            LADA_ASSERT( i_rank_ != i_rank_end_, "Range not iteratable.\n" ) \
            b ## var_iterator result( i_rank_->get<0>().begin() ); \
            for(size_t i(0); i < index_; ++i, ++result); \
            return result; \
          };

      LADA_GETIVAR(const_,const_,const)
      LADA_GETIVAR(,const_,const)
      LADA_GETIVAR(,,const)
      LADA_GETIVAR(,,)
#     undef LADA_GETIVAR


      void CoordRange::rank_range::operator++() 
      {
        if( i_rank_ != i_rank_end_ ) return;
        do 
        { 
          ++i_rank_; 
          if( i_rank_ == i_rank_end_ ) break; 
        } while( i_rank_->get<0>().size() <= index_  ); 
      }
      void CoordRange::const_rank_range::operator++() 
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
