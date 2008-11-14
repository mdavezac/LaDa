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
#include "factory.h"
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
            //! Type of the argument for the general factory.
            typedef TiXmlElement t_XmlNode;
            //! Type of the argument for the attribute factory.
            typedef std::string t_XmlAtt;
            //! Type of the general factories.
            typedef Operators<T_POPULATOR, const t_XmlNode> t_Factory;
            //! Type of the general factories.
            typedef Operators<T_POPULATOR, const t_XmlAtt> t_AttFactory;
            //! Type of the Populator function.
            typedef boost::function<void(T_POPULATOR&)> t_PopFunction;
            //! Type of the this factory.
            typedef XmlOperators<T_INDIVIDUAL, T_POPULATOR> t_This;
          public:
            //! Type of the factory functor.
            typedef boost::function< void( XmlOperators&,
                                     t_PopFunction&, 
                                     const t_XmlNode& ) > t_FactoryFunctor;
            //! Type of the key.
            typedef const std::string t_Key;
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
            void operator()( t_PopFunction& _function, const t_XmlNode &_node );

            //! Sets the default container key.
            void set_default_container_key( const std::string &_functor )
              { default_key_ = _functor; }
           
            //! Connects a functor to the general factory.
            template< class T_FUNCTOR >
              ::LaDa::Factory::ChainConnects<t_This>
                connect( const t_Key& _key, const T_FUNCTOR &_functor );
            //! Connects a functor to the attribute factory.
            template< class T_FUNCTOR >
              void connect_attribute( const t_Key& _key, const T_FUNCTOR &_functor )
              { attfactory_->connect( _key, boost::bind( _functor, *this, _1, _2 ) ); }

            //! Returns true if general factory key exists.
            bool exists( const t_Key& _key ) { return factory_->exists( _key ); }

          protected:
            //! Loops over attributes and performs creations as necessary.
            void try_attribute_factory( t_PopFunction&, const t_XmlNode& );
            //! Returns true if attribute exists.
            bool exists_attribute( const t_Key& _key )
              { return attfactory_->exists( _key ); }

            //! The general factory.
            boost::shared_ptr<t_Factory> factory_;
            //! The attribute factory.
            boost::shared_ptr<t_AttFactory> attfactory_;
            //! A default container creation key.
            std::string default_key_;
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
        void XmlOperators<T_INDIVIDUAL, T_POPULATOR> 
            :: operator()( t_PopFunction& _function, const t_XmlNode &_node )
            { (*this)( default_key_, _function, _node ); }

      template< class T_INDIVIDUAL, class T_POPULATOR > template< class T_FUNCTOR >
        ::LaDa::Factory::ChainConnects< XmlOperators<T_INDIVIDUAL, T_POPULATOR> > 
          XmlOperators<T_INDIVIDUAL, T_POPULATOR>
            :: connect( const t_Key& _key, const T_FUNCTOR &_functor )
            { 
              factory_->connect( _key, boost::bind( _functor, *this, _1, _2 ) ); 
              return ::LaDa::Factory::ChainConnects
                     < 
                       XmlOperators<T_INDIVIDUAL, T_POPULATOR>
                     >( *this );
            }

      template< class T_INDIVIDUAL, class T_POPULATOR > 
        std::ostream& operator<<( std::ostream& _stream, 
                                  const XmlOperators<T_INDIVIDUAL, T_POPULATOR>& _factory )
        {
          _stream << "Xml-tags " << *_factory.factory_
                  << "Xml-attributes " << *_factory.attfactory_;
          return _stream;
        }
    }
  }
}

#endif 
