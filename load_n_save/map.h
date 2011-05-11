//
//  Version: $Id: map.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_MAP_ACTION_H_
#define _LADA_LOADNSAVE_MAP_ACTION_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa 
{
  namespace load_n_save
  {
    //! performs bitwise op on an action.
    template< class T_TYPE >
      class BitwiseOr
      {
        public:
          //! A trait wich says this is an action.
          typedef void action;
          //! The return type is bool.
          typedef bool result_type;
          //! The argument type
          typedef T_TYPE arg_type;
          //! Constructor. 
          BitwiseOr( T_TYPE &_var) : var_(_var) {} 
          //! CopyConstructor. 
          BitwiseOr( BitwiseOr const &_c) : var_(_c.var_ ) {} 
          //! Does bitwise operation.
          result_type operator()( T_TYPE const &_x ) const { var_ |= _x; return true; }

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
