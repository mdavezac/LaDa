//
//  Version: $Id$
//
#ifndef _LADA_OPT_FUSED_FACTORY_H_
#define _LADA_OPT_FUSED_FACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <map>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_front.hpp>

#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/mpl.hpp>

#include <boost/fusion/algorithm/transformation/pop_front.hpp>
#include <boost/fusion/include/pop_front.hpp>
#include <boost/fusion/functional/adapter/unfused_typed.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/fusion/include/invoke.hpp>

#include <boost/function_types/result_type.hpp>
#include <boost/function_types/parameter_types.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/function.hpp>

#include "factory.h"
#include "chainconnects.h"

namespace LaDa
{
  //! Holds factory related objects.
  namespace Factory
  {
    //! \brief Makes pure calls to functors taking no argument and returning no
    //!        values. Functors must be copy constructible.
    template< class T_FUNCTION, class T_KEY >  class FusedFactory;

    //! Dumps a fused factory to a stream.
    template< class T_FUNCTION, class T_KEY >
      std::ostream& operator<<( std::ostream&, const FusedFactory<T_FUNCTION, T_KEY>& );

    template< class T_FUNCTION, class T_KEY = std::string > 
      class FusedFactory 
      {
        friend std::ostream& operator<< <T_FUNCTION, T_KEY>( std::ostream&,
                                                             const FusedFactory<T_FUNCTION, T_KEY>& );
        public:
          //! The type of the key.
          typedef T_KEY t_Key;
          //! Type of the help string.
          typedef std::string t_Help;
          //! Type of the function.
          typedef T_FUNCTION t_Function;
          //! Type of the parameters.
          typedef typename boost::function_types::result_type<t_Function>::type  result_type;
          //! Type of the parameters.
          typedef typename boost::function_types::parameter_types<t_Function>::type 
            t_FuncParameters;

        private:
          //! The type of the map.
          typedef std::map< t_Key, boost::function<t_Function> > t_Map;
          //! Type of the help map.
          typedef std::map< t_Key, t_Help > t_HelpMap;
          //! Type of this class.
          typedef FusedFactory< T_FUNCTION, T_KEY > t_This;
        
        public:
          //! Parameter list with key.
          typedef typename boost::mpl::push_front< t_FuncParameters,
                                                   const t_Key& > :: type t_Parameters;

          //! Constructor.
          FusedFactory() : map_( new t_Map ), helpmap_( new t_HelpMap ) {}
          //! Constructor.
          FusedFactory( const FusedFactory& _c ) : map_( _c.map_ ), helpmap_( _c.helpmap_ ) {}
          //! virtual Destructor.
          virtual ~FusedFactory() {}
 
          //! \brief Adds a new functor. 
          //! \details Throws on duplicates.
          template< class T_FUNCTOR >
            ChainConnects< t_This > connect( const t_Key& _key, const T_FUNCTOR& _functor )
              { return connect<T_FUNCTOR>( _key, "", _functor ); }
          //! \brief Adds a new functor. 
          //! \details Throws on duplicates.
          template< class T_FUNCTOR >
            ChainConnects< t_This > connect( const t_Key& _key, const t_Help& _help, 
                                             const T_FUNCTOR& _functor );

 
          //! performs the call.
          template< class T_SEQUENCE >
            result_type operator()( T_SEQUENCE& );
 
          //! \brief Deletes a connection.
          //! \details Unlike other member functions, this one does not throw if
          //!          \a _key does not exist..
          void disconnect( const t_Key& _key );
          //! Returns true if \a _key exists.
          bool exists( const T_KEY& _key ) const { return map_->end() != map_->find( _key ); }

          //! Returns the help string/object.
          const t_Help& help( const t_Key& _key ) const;
           
        protected:
          //! The map.
          boost::shared_ptr< t_Map > map_;
          //! The map with help stuff.
          boost::shared_ptr< t_HelpMap > helpmap_;
      };

    template< class T_FUNCTION, class T_KEY > template< class T_FUNCTOR >
      ChainConnects< FusedFactory<T_FUNCTION, T_KEY > >
        FusedFactory<T_FUNCTION, T_KEY> :: connect( const t_Key& _key,
                                                    const t_Help& _help,
                                                    const T_FUNCTOR& _functor )
        {
          __DOASSERT( exists( _key ), "Key " << _key << " already exists.\n" )
          typedef typename t_Map::value_type t_type;
          typedef typename t_HelpMap::value_type t_helptype;

          map_->insert( t_type( _key, boost::function< T_FUNCTION >( _functor ) ) );
          helpmap_->insert( t_helptype( _key, _help ) );
          return ChainConnects< t_This >( *this );
        }

    template< class T_FUNCTION, class T_KEY >
      const typename FusedFactory<T_FUNCTION, T_KEY > :: t_Help&
        FusedFactory<T_FUNCTION, T_KEY> :: help( const t_Key& _key ) const
        {
          __DOASSERT( not exists( _key ), "Key " << _key << " does not exists.\n" )
          return (*helpmap_)[_key];
        }


    template< class T_FUNCTION, class T_KEY > template< class T_SEQ >
      inline typename FusedFactory<T_FUNCTION, T_KEY> :: result_type 
        FusedFactory<T_FUNCTION, T_KEY> :: operator()( T_SEQ& _arg )
        {
          const t_Key& key = boost::fusion::at_c<0>( _arg );
          typename t_Map :: iterator i_functor = map_->find( key );
          __DOASSERT( not exists( key ), "Key " << key << " does not exists.\n" )
          boost::fusion::invoke
          ( 
            i_functor->second,
            boost::fusion::pop_front( _arg )
          );
        }

    template< class T_FUNCTION, class T_KEY > 
      inline void FusedFactory<T_FUNCTION, T_KEY> :: disconnect( const t_Key& _key )
      {
         typename t_Map :: iterator i_functor = map_->find( _key );
         if( i_functor != map_->end() ) map_->erase( i_functor );
      }

    template< class T_FUNCTION, class T_KEY >
      std::ostream& operator<<( std::ostream& _stream,
                                const FusedFactory<T_FUNCTION, T_KEY>& _factory )
      {
        typedef typename FusedFactory<T_FUNCTION, T_KEY> :: t_Map :: const_iterator citerator;
        citerator i_map = _factory.map_->begin();
        citerator i_map_end = _factory.map_->end();
        for(; i_map != i_map_end; ++i_map )
          _stream << "  _ \"" << i_map->first << "\" " << _factory.help( i_map->first ) << "\n";
        return _stream;
      }
  }
}

#endif 
