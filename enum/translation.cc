//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "translation.h"
#include "exceptions.h"


namespace LaDa
{

  namespace enumeration
  {
     Translation :: Translation   (math::iVector3d const &_smith, size_t _nsites)
                                : card_(_smith(0)*_smith(1)*_smith(2) * _nsites)
     {
       // loop over possible translations.
       math::iVector3d trans(0,0,0);
       for(trans(0) = 0; trans(0) < _smith(0); ++trans(0))
         for(trans(1) = 0; trans(1) < _smith(1); ++trans(1))
           for(trans(2) = 0; trans(2) < _smith(2); ++trans(2))
           {
             if( trans(0) == 0 and trans(1) == 0 and trans(2) == 0 ) continue;
             std::vector<size_t> permutation; permutation.reserve(card_);
             for(size_t d(0); d<_nsites; ++d)
               for(math::iVector3d x(0, 0,0); x(0) < _smith(0); ++x(0))
                 for(x(1) = 0; x(1) < _smith(1); ++x(1))
                   for(x(2) = 0; x(2) < _smith(2); ++x(2))
                   {
                     math::iVector3d const g
                     (
                       ( trans(0) + x(0) ) % _smith(0),
                       ( trans(1) + x(1) ) % _smith(1),
                       ( trans(2) + x(2) ) % _smith(2) 
                     );
                     permutation.push_back( get_index(d, g, _smith, card_) );
                   }
             permutations_.push_back(permutation);
           }
       i_trans_ = permutations_.rbegin();
       i_trans_end_ = permutations_.rend();
     }

     t_uint Translation::operator()(t_uint _x, FlavorBase const &_flavorbase) const
     {
       if(i_trans_ == i_trans_end_) return _x;
#      ifdef LADA_DEBUG
         if(i_trans_->size() != card_)
           BOOST_THROW_EXCEPTION( internal() << error_string("permutations_ size is incorrect.") );
         if(_flavorbase.size() != card_) BOOST_THROW_EXCEPTION( argument_error());
         if(_x >= _flavorbase.back() * _flavorbase[1]) BOOST_THROW_EXCEPTION( integer_too_large() );
#      endif

       t_uint result(0);
       FlavorBase::const_reverse_iterator i_flavor = _flavorbase.rbegin();
       std::vector<size_t> :: const_iterator i_perm = i_trans_->begin();
       std::vector<size_t> :: const_iterator const i_perm_end = i_trans_->end();
       for(;i_perm != i_perm_end; ++i_flavor, ++i_perm)
       {
         t_uint const flavor( _x / (*i_flavor) );
         _x %= (*i_flavor);

         if(flavor) result += flavor * _flavorbase[*i_perm];
       } // c
       return result;
     }


  } // namespace enumeration.
} // namespace LaDa
