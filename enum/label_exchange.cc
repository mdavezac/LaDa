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
#      ifdef LADA_DEBUG
         if( card_ < 2 ) return _x;
         if( _flavorbase.size() != card_ )
           BOOST_THROW_EXCEPTION( internal() << error_info("_flavorbase size is incorrect.") );
         if( _x >= _flavorbase.back() * _flavorbase[1] )
           BOOST_THROW_EXCEPTION( internal() << error_info("Argument _x is out of range.") );
         if( permutation_.size() != _flavorbase[1] )
           BOOST_THROW_EXCEPTION( internal() << error_info("Incoherent permutation_ or _flavorbase .") );
#      endif
       t_uint result(0);
       FlavorBase::const_reverse_iterator i_flavor = _flavorbase.rbegin();
       for(size_t c(0); c < card_; ++c, ++i_flavor)
       {
         t_uint const flavor( _x / (*i_flavor) );
         _x %= (*i_flavor);

         result += permutation_[flavor] * (*i_flavor);
       } // c
       return result;
     }


  } // namespace enumeration.
} // namespace LaDa
