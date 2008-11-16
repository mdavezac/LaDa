//
//  Version: $Id$
//
#ifndef _LADA_DARWIN_OPERATORS_CALLBACK_H_
#define _LADA_DARWIN_OPERATORS_CALLBACK_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/function.hpp>

//! \cond
template< class EOT > class eoPopulator;
//! \endcond

namespace LaDa
{
  namespace GA
  {
    //! Holds general stuff for GA operators.
    namespace Operator
    {
      //! Creates a unary callback operator.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class UnaryCallBack
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the populator.
            typedef T_POPULATOR t_Populator;
     
            //! Constructor.
            UnaryCallBack() {}
            //! Constructor.
            template< class T_FUNCTOR >
              UnaryCallBack( const T_FUNCTOR& _functor ) : callback_(_functor) {}
            //! Copy Constructor.
            UnaryCallBack( const UnaryCallBack& _c ) : callback_( _c.callback_) {}
            //! Virtual destrcutor, just to make sure.
            virtual ~UnaryCallBack() {}
     
            //! Calls unary with populator.
            void operator()( t_Populator& _populator ) const
            {
              t_Individual& a = *_populator;
              if( callback_( a ) ) a.invalidate();
            }
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator()( t_Individual& _individual ) const
              { return callback_( _individual ); }
            //! Connects the callback.
            template< class T_FUNCTOR > void connect( const T_FUNCTOR& _functor )
              { callback_ = _functor; } 

          private:
            //! The callback object.
            boost::function<bool(t_Individual&) > callback_;
        };

      //! Creates a binary callback operator.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class BinaryCallBack
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the populator.
            typedef T_POPULATOR t_Populator;
     
            //! Constructor.
            BinaryCallBack() {}
            //! Constructor.
            template< class T_FUNCTOR >
              BinaryCallBack( const T_FUNCTOR& _functor ) : callback_(_functor) {}
            //! Copy Constructor.
            BinaryCallBack( const BinaryCallBack& _c ) : callback_( _c.callback_) {}
            //! Virtual destrcutor, just to make sure.
            virtual ~BinaryCallBack() {}
     
            //! Calls unary with populator.
            void operator()( t_Populator& _populator ) const
            {
              t_Individual& a = *_populator;
              const t_Individual& b = _populator.select();
              if( callback_( a, b ) ) a.invalidate();
            }
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator()( t_Individual& _a, const t_Individual& _b ) const
              { return callback_( _a, _b ); }
            //! Connects the callback.
            template< class T_FUNCTOR > void connect( const T_FUNCTOR& _functor )
              { callback_ = _functor; } 
            
          private:
            //! The callback object.
            boost::function<bool(t_Individual&, const t_Individual&) > callback_;
        };

      //! Creates a wrapper around constatnt unary callback operator.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class constUnaryCallBack
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the populator.
            typedef T_POPULATOR t_Populator;
     
            //! Constructor.
            template< class T_FUNCTOR >
              constUnaryCallBack( const T_FUNCTOR& _functor ) : callback_(_functor) {}
            //! Copy Constructor.
            constUnaryCallBack( const constUnaryCallBack& _c ) : callback_( _c.callback_) {}
            //! Virtual destrcutor, just to make sure.
            ~constUnaryCallBack() {}
     
            //! Calls unary with populator.
            void operator()( t_Populator& _populator ) const
            {
              t_Individual& a = *_populator;
              if( callback_( a ) ) a.invalidate();
            }
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator()( t_Individual& _individual ) const
              { return callback_( _individual ); }

          private:
            //! The callback object.
            const boost::function<bool(t_Individual&) > callback_;
        };

      //! Creates a wrapper around a constant binary operator.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class constBinaryCallBack
        {
          public:
            //! Type of the indidividual
            typedef T_INDIVIDUAL t_Individual;
            //! Type of the populator.
            typedef T_POPULATOR t_Populator;
     
            //! Constructor.
            template< class T_FUNCTOR >
              constBinaryCallBack( const T_FUNCTOR& _functor ) : callback_(_functor) {}
            //! Copy Constructor.
            constBinaryCallBack( const constBinaryCallBack& _c ) : callback_( _c.callback_) {}
            //! Virtual destrcutor, just to make sure.
            ~constBinaryCallBack() {}
     
            //! Calls unary with populator.
            void operator()( t_Populator& _populator ) const
            {
              t_Individual& a = *_populator;
              const t_Individual& b = _populator.select();
              if( callback_( a, b ) ) a.invalidate();
            }
            //! Functor over populator. Branches to correct format for t_Derived.
            bool operator()( t_Individual& _a, const t_Individual& _b ) const
              { return callback_( _a, _b ); }
            
          private:
            //! The callback object.
            const boost::function<bool(t_Individual&, const t_Individual&) > callback_;
        };
    }

  }
}



#endif 
