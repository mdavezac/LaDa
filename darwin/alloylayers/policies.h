//
//  Version: $Id$
//
#ifndef _DARWIN_ALLOY_LAYERS_POLICICES_H_
#define _DARWIN_ALLOY_LAYERS_POLICICES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/function.hpp>
#include <boost/signal.hpp>

#include <tinyxml/tinyxml.h>

namespace Crystal { class Structure; }

namespace GA
{
  namespace AlloyLayers
  {
    //! \brief Policy to translate an object containing a vector of reals to an
    //!        epitaxial alloy structure.
    //! \details It is expected that the atoms in the Crystal::Structure
    //!          instances are ordered in direction of growth ( \see
    //!          bool Crystal::create_epitaxial_structure() ). 
    template< class T_OBJECT >
      struct Translate
      {
        //! Type of the Object.
        typedef T_OBJECT t_Object;
        //! From objects to Crystal::Structure. 
        static void translate( const t_Object&, Crystal :: Structure& );
        //! From Crystal::Structure to objects.
        static void translate( const Crystal :: Structure&, t_Object& );
        //! From Crystal::Structure to objects.
        static void translate( const t_Object&, std :: string& );
        //! From Crystal::Structure to objects.
        static void translate( const std::string&, t_Object& );
      };

    //! \brief A function to easily create random individuals.
    //! \details The object of \a _indiv is translated from a random structure.
    template< class T_INDIVIDUAL, class T_TRANSLATE >
      bool initialize( T_INDIVIDUAL &_indiv,
                       Crystal::Structure& _structure,
                       T_TRANSLATE _translate );
                                 
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
            { return (*signal_)( _o, _q ); }
     
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

  }
}

#include "policies.impl.h"

#endif
