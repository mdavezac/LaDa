//
//  Version: $Id$
//
#ifndef _LADA_GA_FACTORY_OPERATORS_H_
#define _LADA_GA_FACTORY_OPERATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
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
      //! \cond
      template< class T_POPULATOR, class T_ARG > class Operators;
      //! \endcond

      //! Dumps keys to stream.
      template< class T_POPULATOR, class T_ARG >
        std::ostream& operator<<( std::ostream&, const Operators<T_POPULATOR, T_ARG>& );

      //! \brief Factory for ga operators.
      //! \details It can take operators of eoPopulators, of one individual,
      //!          and of two individuals.
      template< class T_POPULATOR, class T_ARG = const TiXmlElement >
        class Operators
        {
          friend std::ostream& operator<< <T_POPULATOR, T_ARG>
                                          ( std::ostream&, const Operators<T_POPULATOR, T_ARG>& );
          public:
            //! The type of the key.
            typedef std::string t_Key;
            //! Type of the function.
            typedef boost::function<void(T_POPULATOR& )> t_Function;
          private:
            //! Type of the populator 
            typedef T_POPULATOR t_Populator;
            //! Type of the xml object.
            typedef T_ARG t_Arg;
            //! Pure abstract class for storing purposes.
            class BaseType;
            //! The derived objects which actually store the functors.
            template< class T_FUNCTOR > class DerivedType;
            //! The type of the map.
            typedef boost::ptr_map< t_Key, BaseType > t_Map;
            //! Type of this class.
            typedef Operators< T_POPULATOR, T_ARG > t_This;
          public:
          

            //! Constructor.
            Operators() {}
            //! virtual Destructor.
            virtual ~Operators() {}
      
            //! \brief Adds a new factory function.
            //! \details Throws on duplicate key.
            template< class T_FUNCTOR >
              ::LaDa::Factory::ChainConnects<t_This> connect( const t_Key& _key,
                                                              const T_FUNCTOR& _functor );
            //! performs the call.
            void operator()( const t_Key& _key, t_Function& _function, t_Arg &_arg );
      
            //! \brief Deletes a connection.
            //! \details Unlike other member functions, this one does not throw if
            //!          \a _key does not exist..
            void disconnect( const t_Key& _key );

            //! Returns true if the key exists.
            bool exists( const t_Key& _key )
              { return map_.end() != map_.find( _key ); }
              
             
          protected:
            //! The map.
            t_Map map_;
        };
 
 
    }
  }
}

#include "factory.impl.h"

#endif 
