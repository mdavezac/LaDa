//
//  Version: $Id$
//
#ifndef _LADA_GA_PRINTCALLBACKS_H_
#define _LADA_GA_PRINTCALLBACKS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <boost/function.hpp>
#include <boost/signal.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

#include <tinyxml/tinyxml.h>


namespace LaDa
{

  namespace GA
  {
      //! Policy to assign a printing callback to an object.
      template< class T_OBJECT >
        class PrintCallBack 
        {
          public:
            //! Type of the object.
            typedef T_OBJECT t_Object;
       
            //! Sets the callback.
            template< class T >
            static void connect_print( T _c )
              { callback_ = _c; }

            //! prints using callback.
            static std::ostream& print( std::ostream& _s, const t_Object &_o ) 
              { return _o.callback_( _s, _o ); } 
       
          protected:
            //! The callback object.
            static boost::function<std::ostream&(std::ostream&, const t_Object&) > callback_;
        };
      // A static object to hold a dynamic callback for printing.
      template< class T_OBJECT >
        boost::function< std::ostream& (std::ostream&, const T_OBJECT&) > 
          PrintCallBack<T_OBJECT> :: callback_;
      //! Connects an object to a dynamic print function.
      template< class T_OBJECT >
        std::ostream& operator<<( std::ostream& _s, const PrintCallBack<T_OBJECT>& _o )
        {
          const T_OBJECT& _this = *static_cast<const T_OBJECT*>( &_o );
          return _o.print( _s, _this ); 
        }
      //! Policy to assign a printing callback to an objects.
      template< class T_OBJECT >
        class PrintSignal
        {
          public:
            //! Type of the object.
            typedef T_OBJECT t_Object;
       
            //! Connects a functor/function to the signal.
            static std::ostream& print( std::ostream& _stream, const t_Object& _o )
              { signal_( _stream, _o ); return _stream; }

            //! Connects a functor/function to the signal.
            template< class T_CONNECT >
              static void connect_print( T_CONNECT _c ) { signal_.connect( _c ); }
       
          protected:
            //! Type of the signal.
            typedef boost::signal< void(std::ostream&, const t_Object&)> t_Signal; 
            //! The callback object.
            static t_Signal signal_;
        };
      // A static object to hold a dynamic callback for printing.
      template< class T_OBJECT >
        boost::signal< void (std::ostream&, const T_OBJECT&) > 
          PrintSignal<T_OBJECT> :: signal_;
      //! Connects an object to a dynamic print function.
      template< class T_OBJECT >
        std::ostream& operator<<( std::ostream& _s, const PrintSignal<T_OBJECT>& _o )
        {
          const T_OBJECT& _this = *static_cast<const T_OBJECT*>( &_o );
          return _o.print( _s, _this ); 
        }
  } // namespace GA
} // namespace LaDa

#endif
