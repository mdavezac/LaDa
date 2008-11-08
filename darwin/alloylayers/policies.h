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
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

#include <tinyxml/tinyxml.h>

#include <crystal/layerdepth.h>


namespace LaDa
{
  namespace Crystal { class Structure; }

  namespace GA
  {
    namespace AlloyLayers
    {
      //! A functor which determines whether or not to call an init member function.
      //! \tparam T_POLICY policy for which to check wether is contains an
      //!                  init member function.
      //! \tparam T_CALLER type of the class containing the policy. This type
      //!                  is part of the init member function's signature.
      template< class T_POLICY, class T_CALLER >
        class CallInit
        {
          protected:
            //! \brief A metafunctor using SNIFAE to dectect whether a policy has
            //!        an init method.
            class HasInitMethod
            {
              //! \brief Structure can be instantiated only for certain class.
              //! \details Class U must have a member with the correct signature.
              //!          Function Test then determines the name of this member.
              template<class U, void (U::*)( const T_CALLER& )> struct SFINAE {};
              //! This overload can be instantiated only with the correct class and class member.
              template<class U> static char Test(SFINAE<U, &U::init >*);
              //! This overload can be instantiated with any class.
              template<class U> static int Test(...);
              public:
                //! The resulting value.
                static const bool value = sizeof(Test<T_POLICY>(0)) == sizeof(char);
            };
            //! A simple tag for overloading;
            template< bool docall > class Tag {};

          public:
            //! Type of the policy class.
            typedef T_POLICY t_Policy;
            //! Type of the caller class.
            typedef T_CALLER t_Caller;
            //! \brief Branches through an overloaded function.
            //! \details Casting to policy type is also handled at this point.
            static void call( t_Caller &_caller )
              { branch( _caller, _caller, Tag< HasInitMethod::value >() ); }

          protected:
            //! Does have init member.
            static void branch( t_Policy& _p,
                                const t_Caller & _caller,
                                const Tag< true > & ) { _p.init( _caller ); }
            //! Does not have init member.
            static void branch( t_Policy& _p,
                                const t_Caller & _caller,
                                const Tag< false > & ) { }
        };


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
          //! \brief Initializes the LayerDepth functor.
          //! \note This member is called by the evaluator. It does not need to
          //!       exist in all translate policies, but if it does, it cannot be
          //!       static.
          template< class T_THIS >
            void init( const T_THIS &_this )
             { Translate<T_OBJECT>::layerdepth.set( _this.structure.cell ); }

          protected:
            //! A layer depth operator.
            static Crystal::LayerDepth layerdepth;
        };
      template< class T_OBJECT >
        Crystal::LayerDepth Translate<T_OBJECT> :: layerdepth;

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

      //! \cond 
      namespace details { class isPCallsB {}; }
      namespace details { class isPCallB {}; }
      template< class T_OBJECT > class PrintCallBack;
      //! \endcond
      
      //! Policy to assign a printing callback to an object.
      template< class T_OBJECT >
        class PrintCallBack : public details::isPCallB
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
        typename boost::enable_if
              < 
                boost::is_base_of< details::isPCallB, T_OBJECT >,
                std::ostream
              > :: type & operator<<( std::ostream& _s, const T_OBJECT& _o )
        { return _o.PrintCallBack<T_OBJECT>::print( _s, _o ); }
      //! Policy to assign a printing callback to an objects.
      template< class T_OBJECT >
        class PrintSignal : public details::isPCallsB
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
        typename boost::enable_if
              < 
                boost::is_base_of< details::isPCallsB, T_OBJECT >,
                std::ostream
              > :: type & operator<<( std::ostream& _s, const T_OBJECT& _o )
        { return _o.PrintSignal<T_OBJECT>::print( _s, _o ); }



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

    }
  }
} // namespace LaDa
#include "policies.impl.h"

#endif
