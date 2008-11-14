//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_CALLBACK_H_
#define _LADA_DARWIN_OPERATORS_CALLBACK_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/function.hpp>
#include "discriminate.h"

namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for GA operators.
    namespace Operator
    {
      //! Creates a unary callback operator.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class UnaryCallBack : public PopulatorFunctor
                                     < 
                                       UnaryCallBack<T_INDIVIDUAL, T_POPULATOR>,
                                       T_INDIVIDUAL,
                                       T_POPULATOR
                                     >
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
     
            //! Constructor.
            UnaryCallBack() {}
            //! Copy Constructor.
            UnaryCallBack( const UnaryCallBack& _c ) : callback_( _c.callback_) {}
            //! Virtual destrcutor, just to make sure.
            virtual ~UnaryCallBack() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator()( t_Individual& _individual )
              { return callback_( _individual ); }
            //! Connects the callback.
            template< class T_FUNCTOR > void connect( const T_FUNCTOR& _functor )
              { callback_ = _functor; } 

          private:
            //! The callback object.
            boost::function<bool(t_Individual&) > callback_;
        };

      template< class T_INDIVIDUAL, class T_POPULATOR > 
        struct Discriminate< T_INDIVIDUAL, UnaryCallBack< T_INDIVIDUAL, T_POPULATOR > >
                 : public details::withsizet< Disc::unary_index > {};


      //! Creates a binary callback operator.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class BinaryCallBack : public PopulatorFunctor
                                      < 
                                        BinaryCallBack<T_INDIVIDUAL, T_POPULATOR>,
                                        T_INDIVIDUAL,
                                        T_POPULATOR
                                      >
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
     
            //! Constructor.
            BinaryCallBack() {}
            //! Copy Constructor.
            BinaryCallBack( const BinaryCallBack& _c ) : callback_( _c.callback_) {}
            //! Virtual destrcutor, just to make sure.
            virtual ~BinaryCallBack() {}
     
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator()( t_Individual& _a, const t_Individual& _b )
              { return callback_( _a, _b ); }
            //! Connects the callback.
            template< class T_FUNCTOR > void connect( const T_FUNCTOR& _functor )
              { callback_ = _functor; } 
            
          private:
            //! The callback object.
            boost::function<bool(t_Individual&, const t_Individual&) > callback_;
        };

      template< class T_INDIVIDUAL, class T_POPULATOR > 
        struct Discriminate< T_INDIVIDUAL, BinaryCallBack< T_INDIVIDUAL, T_POPULATOR > >
                 : public details::withsizet< Disc::binary_index > {};

    }

  }
}



#endif 
