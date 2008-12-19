//
//  Version: $Id$
//
#ifndef _LADA_MINIMIZER_DECOUPLED_H_
#define _LADA_MINIMIZER_DECOUPLED_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _DEBUG
#include <stdexcept>
#endif

#include <functional>

#include <tinyxml/tinyxml.h>

#include <opt/debug.h>

#include "any.h"


namespace LaDa
{
  namespace Minimizer
  {
    // \brief Decouples minimization of variables [0, mid[ from [mid, ...[
    class Decoupled 
    {
      public:
        //! Tolerance for each set of variables.
        types::t_real tolerance;
        //! maximum number of iterations
        types::t_int itermax; 
        //! Decouple minimization of variables [0, mid[ from [mid, ...[
        size_t mid; 

        //! Constructor 
        Decoupled() : tolerance(types::tolerance), itermax(500), mid(0) {}
        //! Copy Constructor
        Decoupled   ( const Decoupled &_c ) 
                  : tolerance(_c.tolerance), itermax(_c.itermax), mid(0) {}

        //! Minimization functor
        template< class T_FUNCTION >
          typename T_FUNCTION :: t_Container :: value_type
            operator()( const T_FUNCTION &_func,
                        typename T_FUNCTION :: t_Container &_arg ) const
            { return operator_< T_FUNCTION,
                                typename T_FUNCTION::t_Container,
                                typename T_FUNCTION::t_Container :: value_type
                              >( _func, _arg ); }
        //! Minimization functor
        template< class T_FUNCTION >
          typename T_FUNCTION :: t_Return
            operator()( const T_FUNCTION &_func,
                        typename T_FUNCTION :: t_Arg &_arg ) const
            { return operator_< T_FUNCTION,
                                typename T_FUNCTION :: t_Arg,
                                typename T_FUNCTION :: t_Return
                              >( _func, _arg ); }

        //! Load minimizer parameters from XML
        bool Load( const TiXmlElement &_element );

        //! \brief Loads Minimizer directly from \a _node.
        //! \details If \a _node is not the correct node, the results are undefined.
        bool Load_( const TiXmlElement &_node );

      protected:
        template< class T_FUNCTION, class T_CONTAINER, class T_RETURN >
          T_RETURN operator_( const T_FUNCTION &_func, T_CONTAINER &_arg ) const;

        Any minimizer;
    };

    namespace details
    {
      //! Decoupled functional.
      template< class T_FUNCTION, class T_CONTAINER > class Ranged
      {
        public:
          //! Type of the function;
          typedef T_FUNCTION t_Function;
          //! Type of the argument;
          typedef typename t_Function :: t_Arg t_Arg;
          //! Type of the return;
          typedef typename t_Function :: t_Return t_Return;

          //! Constructor
          Ranged   ( const T_FUNCTION& _func, T_CONTAINER& _cont, int _first, int _last )
                 : function_( _func ), container_( _cont )
          {
            if( _first > 0 ) first_ = size_t( _first );
            if( _first < 0 ) first_ = 0;
            if( _last > 0 ) last_ = size_t( _last );
            if( _last < 0 ) last_ = _cont.size();
          }
          t_Return operator()( t_Return* const _arg ) const
          {
            std::copy( _arg + first_, _arg + last_, container_.begin() + first_ );
            return function_( container_ );
          }
          void gradient( t_Return* const _arg, t_Return *_grad ) const
          {
            std::copy( _arg + first_, _arg + last_, container_.begin() + first_ );
            t_Return grad[ container_.size() ];
            std::fill( grad + first_, grad + last_, t_Return(0) );
            function_.gradient( container_ );
            std::copy( grad + first_, grad + last_, _grad );
          }

        protected:
          //! First in range.
          size_t first_;
          //! last in range.
          size_t last_;
          //! Reference to complete container.
          T_CONTAINER &container_;
          //! Reference to function.
          const T_FUNCTION &function_;
      };
    }
   
    template<typename T_FUNCTION, class T_CONTAINER, class T_RETURN> 
      T_RETURN Decoupled :: operator_( const T_FUNCTION& _function, T_CONTAINER& _arg ) const
      {
        __DOASSERT( mid <= 0, "No need for decoupled minimizer according to input.\n" )
        __DOASSERT( mid >= _arg.size(), "No need for decoupled minimizer according to input.\n" )
        const details::Ranged< T_FUNCTION, T_CONTAINER > funcA( _function, _arg,  -1, mid );
        const details::Ranged< T_FUNCTION, T_CONTAINER > funcB( _function, _arg, mid, -1 );
        T_CONTAINER argA( mid, 0), argB( _arg.size() - mid, 0 );
        std::copy( _arg.begin(), _arg.begin() + mid, argA.begin() );
        std::copy( _arg.begin() + mid, _arg.end(), argB.begin() );
        types::t_real val;
        for( int iter(0); itermax < 0 or iter < itermax; ++iter )
        {
          any( funcA, argA );
          any( funcB, argB );
          if( iter == 0 ) continue;

          types::t_real old = val;
          val = _function( &_arg[0] );
          if( std::abs( old - val ) < tolerance ) break;
        }
        return val;
      }

  }
} // namespace LaDa
#endif
