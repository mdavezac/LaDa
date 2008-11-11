//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_CALLBACK_H_
#define _LADA_DARWIN_OPERATORS_CALLBACK_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "discriminate.h"

namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for GA operators.
    namespace Operator
    {
      //! Creates a unary callback operator.
      template< class T_INDIVIDUAL >
        class UnaryCallBack
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
     
            //! Virtual destrcutor, just to make sure.
            virtual ~UnaryCallBack() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator( t_Individual& _individual )
              { return callback_( _individual ); }
            //! Connects the callback.
            template< class T_FUNCTOR > connect( T_FUNCTOR& _functor )
              { callback_ = _functor; } 
            
          private:
            //! The callback object.
            boost::function<bool(t_Individuals&) > callback_;
        };

      //! Creates a binary callback operator.
      template< class T_INDIVIDUAL >
        class BinaryCallBack
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
     
            //! Virtual destrcutor, just to make sure.
            virtual ~UnaryCallBack() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator( t_Individual& _a, const t_Individual& _b )
              { return callback_( _a, _b ); }
            //! Connects the callback.
            template< class T_FUNCTOR > connect( T_FUNCTOR& _functor )
              { callback_ = _functor; } 
            
          private:
            //! The callback object.
            boost::function<bool(t_Individuals&, const t_Individual&) > callback_;
        };


    }

  }
}



#endif 
