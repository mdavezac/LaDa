//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "label_exchange.h"
#include "exceptions.h"


namespace LaDa
{

  namespace enumeration
  {
     t_uint LabelExchange::operator()(t_uint _x, FlavorBase const &_flavorbase) const
     {
       if( card_ < 2 ) return _x;
#      ifdef LADA_DEBUG
         if(_flavorbase.size() != card_)
           BOOST_THROW_EXCEPTION( internal() << error_string("_flavorbase size is incorrect.") );
         if(_flavorbase.size() != card_) BOOST_THROW_EXCEPTION( argument_error());
         if(_x >= _flavorbase.back() * _flavorbase[1]) BOOST_THROW_EXCEPTION( integer_too_large() );
#      endif
       t_uint result(0);
       FlavorBase::const_reverse_iterator i_flavor = _flavorbase.rbegin();
       for(size_t c(0); c < card_; ++c, ++i_flavor)
       {
         t_uint const flavor( permutation_[_x / (*i_flavor)] );
         _x %= (*i_flavor);

         if( flavor ) result += flavor * (*i_flavor);
       } // c
       return result;
     }


  } // namespace enumeration.
} // namespace LaDa
