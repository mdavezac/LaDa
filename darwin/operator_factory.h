//
//  Version: $Id$
//
#ifndef _LADA_OPT_FACTORY_H_
#define _LADA_OPT_FACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <opt/factory.h>

namespace LaDa
{
  namespace GA
  {
    //! Holds  GA factory related objects.
    namespace Factory
    {
      //! \brief Factory for ga operators.
      //! \details It can take operators of eoPopulators, of one individual,
      //!          and of two individuals.
      template< class T_INDIVIDUAL >
        class Operators : public opt::ChainConnects< T_INDIVIDUAL >
        {
          //! Pure abstract class for storing purposes.
          class BaseType;
          //! The derived objects which actually store the functors.
          template< class T_FUNCTOR > class DerivedType;
          //! The type of the key.
          typedef const std::string t_Key;
          //! The type of the map.
          typedef boost::ptr_map< t_Key, BaseType > t_Map;
          //! Wraps an operator over an individual.
          template< class T_FUNCTOR > class WrapMonOp;
          //! Wraps an operator over two individuals.
          template< class T_FUNCTOR > class WrapBinOp;
          //! Discriminates between operator types.
          class IsMonary;
          //! Discriminates between operator types.
          class IsGenOp;
          
          public:
            //! Constructor.
            GAOperators( eoState _eostates ) : eostates_( _eostates ){}
            //! virtual Destructor.
            virtual ~GAOperators() {}
      
            //! \brief Adds a new factory function.
            //! \details Throws on duplicate key.
            template< class T_FUNCTOR >
              ChainConnects<T_INDIVIDUAL> connect( const t_Key& _key, T_FUNCTOR* _functor );
            //! \brief Adds a new factory function.
            //! \details Throws on duplicate key.
            template< class T_FUNCTOR >
      
            //! performs the call.
            eoGenOp<t_Individual>* operator()( const t_Key& _key )
      
            //! \brief Deletes a connection.
            //! \details Unlike other member functions, this one does not throw if
            //!          \a _key does not exist..
            void disconnect( const t_Key& _key );
             
          protected:
            //! The map.
            t_Map map_;
            //! The garbage collector.
            eoState& eostates_;
        };
 
      template< class T_INDIVIDUAL >
        struct PureCalls<T_INDIVIDUAL> :: BaseType 
        {
          //! The virtual functor. 
          virtual void create_basetype( eoGenOp<t_Individual >* ) = 0;
          //! A virtual Destructor.
          virtual ~BaseType() {};
        };
      template< class T_FUNCTOR >
        class PureCalls :: DerivedType : public BaseType 
        {
          public:
            //! the type of the functor.
            typedef T_FUNCTOR t_Functor;
            //! Constructor
            DerivedType   ( const t_Functor& _functor )
                        : functor_( new t_Functor( _functor ) ) {}
            //! Copy Constructor.
            DerivedType( const DerivedType& _c ) : functor_( _c.functor_ ) {}
            //! Destructor.
            virtual ~DerivedType() {}
           
            //! The virtual functor. 
            void create_basetype( eoGenOpt<t_Individual>* _result ) { (*functor)( _result ); }
 
          protected:
            //! Holds the functor.
            boost::shared_ptr< t_Functor > functor_;
        };
 
      template< class T_FUNCTOR >
        Operators<T_INDIVIDUAL> :: ChainConnects 
          Operators<T_INDIVIDUAL> :: connect( const t_Key& _key,
                                              const T_FUNCTOR& _functor )
          {
            __DOASSERT( map_.end() != map_.find( _key ),
                        "Key " << _key << " already exists.\n" )
            typedef typename t_Map :: value_type value_type;
            map_.insert( _key, new DerivedType<T_FUNCTOR>( _functor ) );
            return ChainConnects( *this );
          }

      template< class T_INDIVIDUAL >
        eoGenOp<T_INDIVIDUAL>* Operators<T_INDIVIDUAL> :: operator()( const t_Key& _key )
        {
           t_Map :: iterator i_functor = map_.find( _key );
           __DOASSERT( i_functor == map_.end() , "Key " << _key << " does not exists.\n" )
           eoGenOp<t_Individual>* result;
           (*i_functor->second)( result, *this );
           eostates_.storeFunctor( result );
           return result;
        }
 
      template< class T_INDIVIDUAL >
        void Operators<T_INDIVIDUAL> :: disconnect( const t_Key& _key )
        {
           t_Map :: iterator i_functor = map_.find( _key );
           if( i_functor != map_.end() ) map_.erase( i_functor );
        }

 
    }
  }
}

#include "operator_factory.impl.h"

#endif 
