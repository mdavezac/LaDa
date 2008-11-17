//
//  Version: $Id$
//
#ifndef _LADA_GA_FACTORY_XMLOPERATORS_H_
#define _LADA_GA_FACTORY_XMLOPERATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <opt/factory.h>
#include <opt/chainconnects.h>

#include "container.h"

namespace LaDa
{
  namespace GA
  {
    //! Holds  GA factory related objects.
    namespace Factory
    {
      //! \cond
      template< class T_INDIVIDUAL, class T_POPULATOR > class XmlOperators;
      //! \endcond
 
      //! Dumps keys to stream.
      template< class T_INDIVIDUAL, class T_POPULATOR >
        std::ostream& operator<<( std::ostream&, const XmlOperators<T_INDIVIDUAL, T_POPULATOR>& );

      //! Reads from xml and creates operators.
      template< class T_INDIVIDUAL, class T_POPULATOR = eoPopulator<T_INDIVIDUAL> >
        class XmlOperators
        {
          friend std::ostream& operator<< <T_INDIVIDUAL, T_POPULATOR>
                                          ( std::ostream&,
                                            const XmlOperators<T_INDIVIDUAL, T_POPULATOR>& );
          public:
            //! Type of the general key.
            typedef std::string t_Key;
            //! Type of the attribute key.
            typedef std::string t_attKey;
          protected:
            //! Type of the argument for the general factory.
            typedef TiXmlElement t_XmlNode;
            //! Type of the argument for the attribute factory.
            typedef std::string t_XmlAtt;
            //! Type of the Populator function.
            typedef boost::function<void(T_POPULATOR&)> t_PopFunction;
            //! Type of the general factories.
            typedef ::LaDa::Factory::Factory< void(t_PopFunction&, const t_XmlNode&),
                                              t_Key > t_Factory;
            //! Type of the general factories.
            typedef ::LaDa::Factory::Factory< void(t_PopFunction&, const t_XmlAtt&),
                                              t_attKey > t_AttFactory;
            
            //! Type of the this factory.
            typedef XmlOperators<T_INDIVIDUAL, T_POPULATOR> t_This;

            //! Type of the connect return.
            typedef ::LaDa::Factory::ChainConnects< t_This > t_ConnectReturn;
            //! Type of the connect return.
            typedef void t_AttConnectReturn;
          public:
            //! Type of the help string.
            typedef typename t_Factory :: t_Help t_Help;
            //! Type of the factory functor.
            typedef boost::function< void( XmlOperators&,
                                     t_PopFunction&, 
                                     const t_XmlNode& ) > t_FactoryFunctor;
            //! Type of the populator.
            typedef T_POPULATOR t_Populator;
            //! Type of the individual.
            typedef T_INDIVIDUAL t_Individual;

            //! Constructor.
            XmlOperators() : factory_( new t_Factory ), attfactory_( new t_AttFactory ),
                             default_key_( "And" ) {}
            //! Copy Constructor.
            XmlOperators   ( const XmlOperators& _c ) 
                         : factory_( _c.factory_ ), attfactory_( _c.attfactory_ ),
                           default_key_( _c.default_key_ ) {}

            //! Searches amongst general factory for a function key an creates it.
            void operator()( const t_Key &, t_PopFunction&, const t_XmlNode & );
            //! Creates an operator from current node using the default container key.
            void operator()( t_PopFunction& _function, const t_XmlNode &_node )
              { (*this)( default_key_, _function, _node ); }

            //! Sets the default container key.
            void set_default_container_key( const t_Key &_functor )
              { default_key_ = _functor; }
           
            //! Connects a functor to the general factory.
            template< class T_FUNCTOR >
              t_ConnectReturn connect( const t_Key& _key, const T_FUNCTOR &_functor )
                { return connect( _key, "", _functor ); }
            //! Connects a functor to the general factory.
            template< class T_FUNCTOR >
              t_ConnectReturn connect( const t_Key& _key, const t_Help& _help,
                                       const T_FUNCTOR &_functor )
              {
                factory_->connect( _key, _help, boost::bind( _functor, *this, _1, _2 ) ); 
                return t_ConnectReturn( *this );
              }
            //! Connects a functor to the attribute factory.
            template< class T_FUNCTOR >
              t_AttConnectReturn connect_attribute( const t_attKey & _key,
                                                    const T_FUNCTOR &_functor )
              { connect_attribute( _key, " ", boost::bind( _functor, *this, _1, _2 ) ); }
            template< class T_FUNCTOR >
              t_AttConnectReturn connect_attribute( const t_attKey & _key, const t_Help& _help,
                                                    const T_FUNCTOR &_functor )
              { attfactory_->connect( _key, _help, boost::bind( _functor, *this, _1, _2 ) ); }

            //! Returns true if general factory key exists.
            bool exists( const t_Key& _key ) { return factory_->exists( _key ); }

          protected:
            //! Loops over attributes and performs creations as necessary.
            void try_attribute_factory( t_PopFunction&, const t_XmlNode& );
            //! Returns true if attribute exists.
            bool exists_attribute( const t_attKey& _key )
              { return attfactory_->exists( _key ); }

            //! The general factory.
            boost::shared_ptr<t_Factory> factory_;
            //! The attribute factory.
            boost::shared_ptr<t_AttFactory> attfactory_;
            //! A default container creation key.
            t_Key default_key_;
        };

      template< class T_INDIVIDUAL, class T_POPULATOR >
        void XmlOperators<T_INDIVIDUAL, T_POPULATOR> 
          :: try_attribute_factory( t_PopFunction& _function, const t_XmlNode& _node )
          {
            const TiXmlAttribute* att = _node.FirstAttribute(); 
            for(; att; att = att->Next() )
            {
              const std::string key( att->Name() );
              if( not exists_attribute( key ) ) continue;
              const std::string value( att->Value() );
              (*attfactory_)( key, _function, value );
            }
          }
      template< class T_INDIVIDUAL, class T_POPULATOR >
        void XmlOperators<T_INDIVIDUAL, T_POPULATOR> 
            :: operator()( const t_Key &_key, t_PopFunction& _function, const t_XmlNode & _node )
            {
              __DOASSERT( not exists( _key ), _key << " is not a registered operator.\n" )
              (*factory_)( _key, _function, _node );
              try_attribute_factory( _function, _node );
            }


      template< class T_INDIVIDUAL, class T_POPULATOR > 
        std::ostream& operator<<( std::ostream& _stream, 
                                  const XmlOperators<T_INDIVIDUAL, T_POPULATOR>& _factory )
        {
          _stream << "Xml-tags: \n" << *_factory.factory_
                  << "Xml-attributes: \n" << *_factory.attfactory_;
          return _stream;
        }
    }
  }
}

#endif 
