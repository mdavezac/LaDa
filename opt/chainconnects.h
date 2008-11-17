//
//  Version: $Id$
//
#ifndef _LADA_OPT_FACTORY_CHAINCONNECTS_H_
#define _LADA_OPT_FACTORY_CHAINCONNECTS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace Factory
  {
    //! Helps chaining together calls to T_THIS :: connect().
    template< class T_THIS > 
      class ChainConnects
      {
        public:
          //! Type of the derived class.
          typedef T_THIS t_This;
         
     
          //! Copy constructor. 
          ChainConnects( const ChainConnects &_c ) : this_( _c.this_ ) {}
          //! Constructor. 
          ChainConnects( t_This &_this ) : this_( _this ) {}

          //! Functor which chains calls to connect.
          template< class T_FUNCTOR >
            ChainConnects& operator()( const typename t_This::t_Key& _key,
                                       const T_FUNCTOR& _functor )
              { this_.connect( _key, _functor ); return *this; }
          //! Functor which chains calls to connect.
          template< class T_FUNCTOR >
            ChainConnects& operator()( const typename t_This::t_Key& _key,
                                       const typename t_This::t_Help& _help,
                                       const T_FUNCTOR& _functor )
              { this_.connect( _key, _help, _functor ); return *this; }
        private:
          //! Holds a reference to PureCalls.
          t_This &this_;
      };
   }
}

#endif

