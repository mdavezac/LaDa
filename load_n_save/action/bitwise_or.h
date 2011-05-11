//
//  Version: $Id: bitwise_or.h 1250 2009-07-26 21:04:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_BITWISE_OR_ACTION_H_
#define _LADA_LOADNSAVE_BITWISE_OR_ACTION_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "initializer/parse_value.h"

namespace LaDa 
{
  namespace load_n_save
  {
    //! performs bitwise or on an action.
    template< class T_TYPE >
      class BitwiseOr
      {
        public:
          //! A trait wich says this is an action.
          typedef void action;
          //! The return type is bool.
          typedef bool result_type;
          //! The argument type
          typedef std::string arg_type;
          //! Constructor. 
          BitwiseOr( T_TYPE &_var) : var_(_var) {} 
          //! CopyConstructor. 
          BitwiseOr( BitwiseOr const &_c) : var_(_c.var_ ) {} 
          //! Does bitwise operation.
          result_type operator()( std::string const &_x ) const
          {
            T_TYPE arg;
            if( not initializer::parse_value( _x, arg ) ) return false;
            var_ |= arg; 
            return true; 
          }

        protected:
          //! Holds reference to variable.
          T_TYPE &var_;
      };

    //! Returns a bitwise or action.
    template< class T_TYPE >
      BitwiseOr<T_TYPE> bitwise_or( T_TYPE &_a ) { return BitwiseOr<T_TYPE>(_a); }

  }
}

#endif
