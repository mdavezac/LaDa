//
//  Version: $Id$
//
#ifndef _LADA_FACTORY_FACTORY_H_
#define _LADA_FACTORY_FACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <map>

#include "fusedfactory.h"

namespace LaDa
{
  //! Holds factory related objects.
  namespace Factory
  {
    //! The callable factory.
    template< class T_FUNCTION, class T_KEY >  class Factory;

    //! Dumps a factory to a stream.
    template< class T_FUNCTION, class T_KEY >
      std::ostream& operator<<( std::ostream& _stream,
                                const Factory<T_FUNCTION, T_KEY>& _factory );
    namespace details
    {
      //! An intermediate base class for correct constructor call ordering.
      template< class T_FUNCTION, class T_KEY >
        class Intermediate  
        {
            //! Type of the fused class.
            typedef FusedFactory<T_FUNCTION, T_KEY> t_Fused;
            //! Type of the base class.
            typedef boost::fusion::unfused_typed< t_Fused, typename t_Fused :: t_Parameters> t_Base; 
          public:
            //! Constructor.
            Intermediate() {}
            //! Copy Constructor
            Intermediate( const Intermediate& _c ) : fused_(_c.fused_) {}
      
            //! Forwards to FusedFactory :: connect()
            template< class T_FUNCTOR >
              ChainConnects< t_Fused > connect( const T_KEY& _key, const T_FUNCTOR& _functor )
                { return fused_. connect<T_FUNCTOR>( _key, _functor ); }
            //! Forwards to FusedFactory :: connect()
            template< class T_FUNCTOR >
              ChainConnects< t_Fused > connect( const T_KEY& _key,
                                                const typename t_Fused::t_Help& _help,
                                                const T_FUNCTOR& _functor )
                { return fused_. connect<T_FUNCTOR>( _key, _help, _functor ); }
      
            //! Forwards to FusedFactory :: disconnect()
            void disconnect( const T_KEY& _key )
              { fused_.disconnect( _key ); }
            //! Returns true if \a _key exists.
            bool exists( const T_KEY& _key ) const { return fused_.exists( _key ); }
            //! Returns help string.
            const typename t_Fused::t_Help& help( const T_KEY& _key ) const
              { return fused_.help( _key ); }
      
            protected:
             t_Fused fused_;
        };
    }
   
    //! Creates a factory with an unfused operator().
    template< class T_FUNCTION, class T_KEY = std::string >
      class Factory : public details :: Intermediate< T_FUNCTION, T_KEY >,
                      public boost::fusion::unfused_typed
                              <
                                const FusedFactory<T_FUNCTION, T_KEY>,
                                typename FusedFactory<T_FUNCTION, T_KEY> :: t_Parameters 
                              > 
      {
          friend std::ostream& operator<< <T_FUNCTION, T_KEY>
                                       ( std::ostream& _stream,
                                         const Factory<T_FUNCTION, T_KEY>& _factory );
          //! Type of the fused class.
          typedef FusedFactory<T_FUNCTION, T_KEY> t_Fused;
          //! Type of the fused base class.
          typedef details::Intermediate< T_FUNCTION, T_KEY > t_Intermediate;
          //! Type of the unfused base class.
          typedef boost::fusion::unfused_typed< const t_Fused,
                                                typename t_Fused :: t_Parameters> t_UnFusedBase; 
        public:
          //! The function type.
          typedef T_FUNCTION t_Function;
          //! The key type.
          typedef T_KEY t_Key;
          //! The help type.
          typedef typename t_Fused :: t_Help t_Help;
          //! Type of the connect return.
          typedef ChainConnects< t_Fused > t_ConnectReturn;
          //! Constructor.
          Factory() : t_Intermediate(), t_UnFusedBase( t_Intermediate::fused_ ) {}
          //! Copy Constructor.
          Factory   ( const Factory& _c )
                  : t_Intermediate( _c ), t_UnFusedBase( t_Intermediate::fused_ ) {}
      };


    template< class T_FUNCTION, class T_KEY >
      std::ostream& operator<<( std::ostream& _stream,
                                const Factory<T_FUNCTION, T_KEY>& _factory )
      { return _stream << _factory.details::Intermediate<T_FUNCTION,T_KEY>::fused_; }

  }
}

#endif 
