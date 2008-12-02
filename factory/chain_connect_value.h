//
//  Version: $Id: chainconnects.h 860 2008-11-17 18:37:10Z davezac $
//
#ifndef _LADA_FACTORY_CHAINCONNECTVALUE_H_
#define _LADA_FACTORY_CHAINCONNECTVALUE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace Factory
  {
    //! Helps chaining together calls to T_THIS :: connect().
    template< class T_THIS > 
      class ChainConnectValue
      {
        public:
          //! Type of the derived class.
          typedef T_THIS t_This;
         
     
          //! Copy constructor. 
          ChainConnectValue( const ChainConnectValue &_c ) : this_( _c.this_ ) {}
          //! Constructor. 
          ChainConnectValue( t_This &_this ) : this_( _this ) {}

          //! Functor which chains calls to connect.
          template< class T_FUNCTOR >
            ChainConnectValue& operator()( const typename t_This::t_Key& _key,
                                           const T_FUNCTOR& _functor )
              { this_.connect_value( _key, _functor ); return *this; }
          //! Functor which chains calls to connect.
          template< class T_FUNCTOR >
            ChainConnectValue& operator()( const typename t_This::t_Key& _key,
                                           const typename t_This::t_Help& _help,
                                           const T_FUNCTOR& _functor )
              { this_.connect_value( _key, _help, _functor ); return *this; }
        private:
          //! Holds a reference to PureCalls.
          t_This &this_;
      };
   }
}

#endif

