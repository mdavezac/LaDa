//
//  Version: $Id$
//
#ifndef _LADA_GA_ASSIGNCALLBACKS_H_
#define _LADA_GA_ASSIGNCALLBACKS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/function.hpp>
#include <boost/signal.hpp>
#include <boost/utility/enable_if.hpp>

#include <tinyxml/tinyxml.h>

namespace LaDa
{
  namespace GA
  {
      //! \brief Policy to connect a single callback function assigning
      //!        values from an object to a GA quantity/raw fitness.
      template< class T_OBJECT, class T_QUANTITIES >
        class AssignCallBack
        {
          public:
            //! Type of the object.
            typedef T_OBJECT t_Object;
            //! Type of the quantities.
            typedef T_QUANTITIES t_Quantities;
       
            //! Assigns a value from object to a quantity.
            void assign( const t_Object& _o, t_Quantities &_q ) const
              { return callback_( _o, _q ); }
       
            //! Sets the callback.
            void connect( void (*_c)( t_Object, t_Quantities ) )
              { callback_ = _c; }
       
          protected:
            //! The callback object.
            boost::function<void(const t_Object&, t_Quantities&) > callback_;
        };


      //! \brief Policy to connect any number of callback functions assigning
      //!        values from an object to a GA quantity/raw fitness.
      template< class T_OBJECT, class T_QUANTITIES >
        class AssignSignal
        {
          public:
            //! Type of the object.
            typedef T_OBJECT t_Object;
            //! Type of the quantities.
            typedef T_QUANTITIES t_Quantities;

            //! Constructor.
            AssignSignal() : signal_( new t_Signal ) {}
            //! Copy Constructor.
            AssignSignal( const AssignSignal & _c ) : signal_( _c.signal_ ) {}
       
            //! Connects a functor/function to the signal.
            void assign( const t_Object& _o, t_Quantities &_q ) const
              { _q.clear(); return (*signal_)( _o, _q ); }
       
            //! Connects a functor/function to the signal.
            template< class T_CONNECT >
              void connect( T_CONNECT _c ) { signal_->connect( _c ); }
       
          protected:
            //! Type of the signal.
            typedef boost::signal< void(const t_Object&, t_Quantities&)> t_Signal; 
            //! \brief The signal callback container.
            //! \details This is a boost shared pointer, and as such, signal_ is
            //!          copyable.
            boost::shared_ptr< t_Signal > signal_;
        };

  } // namespace GA
} // namespace LaDa

#endif
