//
//  Version: $Id: xmlfactory.h 860 2008-11-17 18:37:10Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_FACTORY_H_
#define _LADA_GA_CHECKPOINTS_FACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/static_assert.hpp>

#include <factory/factory.h>
#include <factory/chainconnects.h>

#include <factory/xmlkey.h>

namespace LaDa
{
  namespace GA
  {
    //! Holds  GA factory related objects.
    namespace Factory
    {
      //! \cond
      template< class T_INDIVIDUAL, class T_POPULATOR > class Checkpoints;
      //! \endcond
 
      //! Dumps keys to stream.
      template< class T_INDIVIDUAL, class T_POPULATOR >
        std::ostream& operator<<( std::ostream&, const CheckPoints<T_INDIVIDUAL, T_POPULATOR>& );

      //! Reads from xml and creates operators.
      template< class T_CHECKPOINT, class T_KEY = XmlKey >
        class CheckPoints
        {
            friend std::ostream& operator<< <T_INDIVIDUAL, T_POPULATOR>
                                            ( std::ostream&,
                                              const CheckPoints<T_INDIVIDUAL, T_POPULATOR>& );
            //! Type of the this factory.
            typedef CheckPoints<T_INDIVIDUAL, T_POPULATOR> t_This;
            //! Type of the argument for the general factory.
            typedef TiXmlElement t_XmlNode;
            //! Type of the argument for the attribute factory.
            typedef std::string t_XmlAtt;
            //! Type of the connect return.
            typedef ::LaDa::Factory::ChainConnects< t_This > t_ConnectReturn;
            //! Type of the connect return.
            typedef void t_AttConnectReturn;
          public:
            //! Type of the checkpoint class.
            typedef T_CHECKPOINT t_CheckPoint;
            //! Type of the general key.
            typedef std::string t_Key;
            //! Type of node function.
            typedef void( t_CheckPoint&, const t_XmlNode& ) t_NodeFunction;
            //! Type of attribute function.
            typedef void( t_CheckPoint&, const t_XmlAttribute& ) t_AttributeFunction;
            //! Type of value function.
            typedef void( t_CheckPoint& ) t_ValueFunction;

          public:
            //! Type of the help string.
            typedef typename t_Factory :: t_Help t_Help;
            //! Type of the factory functor.
            typedef boost::function< void( CheckPoints&,
                                     t_PopFunction&, 
                                     const t_XmlNode& ) > t_FactoryFunctor;

            //! Constructor.
            CheckPoints() : node_factory_( new t_NodeFactory ),
                            attribute_factory_( new t_AttributeFactory ),
                            value_factory_( new t_ValueFactory ) {}
            //! Copy Constructor.
            CheckPoints   ( const CheckPoints& _c ) 
                        : node_factory_( _c.node_factory_ ),
                          attribute_factory_( _c.attribute_factory_ ),
                          value_factory_( _c.value_factory_ ) {}

            //! Searches amongst all factories for a function key and creates it.
            void operator()( t_CheckPoint&, const t_XmlNode & );

            //! Connects a functor to the general factory.
            template< class T_FUNCTOR >
              t_ConnectReturn connect( const t_Key& _key, const T_FUNCTOR &_functor )
                { return connect( _key, "", _functor ); }
            //! Connects a functor to the general factory.
            template< class T_FUNCTOR >
              t_ConnectReturn connect( const t_Key& _key, const t_Help& _help,
                                       const T_FUNCTOR &_functor );
            //! Returns true if factory key exists.
            bool exists( const t_Key& _key )
            { return     node_factory_->exists( _key ) 
                     or attribute_factory_->exists( _key )
                     or value_factory_->exists( _key ); }

          protected:
            //! Helper class to check conversion to boost::function<t_Continuator>
            template< class T_FUNCTOR >
              struct is_node : public boost::is_convertible
                        < boost::function< t_NodeFunction >, T_FUNCTOR > {};
            //! Helper class to check conversion to boost::function<t_Updater>
            template< class T_FUNCTOR >
              struct is_attribute : public boost::is_convertible
                        < boost::function< t_NodeFunction >, T_FUNCTOR > {};
            //! Helper class to check conversion to boost::function<t_Statistics>
            template< class T_FUNCTOR >
              struct is_attribute : public boost::is_convertible
                        < boost::function< t_Attribute >, T_FUNCTOR > {};

            //! connects an node function.
            template< class T_FUNCTOR >
              typename boost::enable_if< is_node< T_FUNCTOR > > :: type
                connect_( const t_Key& _key, const t_Help& _help, const T_FUNCTOR& _functor )
                  { node_factory_.connect( _key, _help, _functor ); }
            //! connects an attribute function.
            template< class T_FUNCTOR >
              typename boost::enable_if< is_attribute< T_FUNCTOR > > :: type
                connect_( const t_Key& _key, const t_Help& _help, const T_FUNCTOR& _functor )
                  { attribute_factory_.connect( _key, _help, _functor ); }
            //! connects a value function.
            template< class T_FUNCTOR >
              typename boost::enable_if< is_value< T_FUNCTOR > > :: type
                connect_( const t_Key& _key, const t_Help& _help, const T_FUNCTOR& _functor )
                  { value_factory_.connect( _key, _help, _functor ); }

            //! Type of the node factory.
            typedef ::LaDa::Factory::Factory< t_NodeFunction, t_Key > t_NodeFactory;
            //! Type of the node factory.
            typedef ::LaDa::Factory::Factory< t_AttributeFunction, t_Key > t_AttributeFactory;
            //! Type of the node factory.
            typedef ::LaDa::Factory::Factory< t_ValueFunction, t_Key > t_ValueFactory;

            //! The node factory.
            boost::shared_ptr<t_Factory> node_factory_;
            //! The attribute factory.
            boost::shared_ptr<t_AttFactory> attribute_factory_;
            //! The value factory.
            boost::shared_ptr<t_AttFactory> value_factory_;
        };

      template< class T_CHECKPOINT, class T_KEY > 
        std::ostream& operator<<( std::ostream& _stream, 
                                  const CheckPoints<T_CHECKPOINT, T_KEY>& _factory )
        {
          _stream << *_factory.node_factory_ << "\n"
                  << *_factory.attribute_factory_ << "\n"
                  << *_factory.value_factory_;
          return _stream;
        }

      template< class T_CHECKPOINT, class T_KEY > template <class T_FUNCTOR>
        typename CheckPoints<T_CHECKPOINT, T_KEY> :: t_ConnectReturn
          CheckPoints<T_CHECKPOINT, T_KEY > :: connect( const t_Key& _key, 
                                                        const t_Help& _help,
                                                        const T_FUNCTOR& _functor )
          {
            BOOST_STATIC_ASSERT(    is_node<T_FUNCTOR>::value 
                                 or is_attribute<T_FUNCTOR>::value 
                                 or is_value<T_FUNCTOR> :: value )
            connect_( _key, _help, _functor );
          }

    }
  }
}

#endif 
