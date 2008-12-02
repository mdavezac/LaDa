//
//  Version: $Id: xmlfactory.h 860 2008-11-17 18:37:10Z davezac $
//
#ifndef _LADA_FACTORY_XMLFUSEDFACTORY_H_
#define _LADA_FACTORY_XMLFUSEDFACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/equal.hpp>

#include "fusedfactory.h"
#include "chainconnects.h"
#include "functiontype.h"
#include "chain_connect_node.h"
#include "chain_connect_attribute.h"
#include "chain_connect_value.h"

#include "xmlkey.h"

namespace LaDa
{
  //! Holds  GA factory related objects.
  namespace Factory
  {
    //! \cond
    template< class T_FUNCTION, class T_KEY > class XmlFusedFactory;
    //! \endcond
 
    //! Dumps keys to stream.
    template< class T_FUNCTION, class T_KEY >
      std::ostream& operator<<( std::ostream&, const XmlFusedFactory<T_FUNCTION, T_KEY>& );

    namespace details
    {
      //! A class to disable implicit conversion from std::string to XmlElement.
      class XmlElement : public TiXmlElement
      {
        public:
          //! Constructor.
          XmlElement( const TiXmlElement& _element ) : TiXmlElement( _element ) {}
          //! Copy Constructor.
          XmlElement( const XmlElement& _element ) : TiXmlElement( _element ) {}
          //! Destructor.
          ~XmlElement() {}
        private:
          //! Private copy constructor.
          XmlElement( const std::string& ) : TiXmlElement("") {}
      };
    }

    //! Reads from xml and creates operators.
    template< class T_FUNCTION, class T_KEY = XmlKey >
      class XmlFusedFactory
      {
          friend std::ostream& operator<< <T_FUNCTION, T_KEY>
                                          ( std::ostream&,
                                            const XmlFusedFactory<T_FUNCTION, T_KEY>& );
          //! Type of the this factory.
          typedef XmlFusedFactory<T_FUNCTION, T_KEY> t_This;
        protected:
          //! Type of the argument for the general factory.
          typedef details::XmlElement t_XmlElement;
          //! Type of the argument for the attribute factory.
          typedef std::string t_XmlAttribute;
          //! Type returned by connect_node().
          typedef ChainConnectNode<t_This> t_NodeReturn;
          //! Type returned by connect_attribute().
          typedef ChainConnectAttribute<t_This> t_AttributeReturn;
          //! Type returned by connect_attribute().
          typedef ChainConnectValue<t_This> t_ValueReturn;
        public:
          //! Type of the functions.
          typedef T_FUNCTION t_Function;
          //! Type of the general key.
          typedef T_KEY t_Key;
          //! \brief Type of the help string.
          //! \todo bind this type directly to Factory::Factory.
          typedef std::string t_Help;
          //! Type of the return.
          typedef typename boost::function_types::result_type<t_Function>::type  result_type;
          //! Type of node function.
          typedef typename push_back_arg< t_Function, const t_XmlElement& > :: type t_NodeFunction;
          //! Type of attribute function.
          typedef typename push_back_arg< t_Function, const t_XmlAttribute& > :: type
            t_AttributeFunction;
          //! Type of value function.
          typedef T_FUNCTION t_ValueFunction;

          BOOST_STATIC_ASSERT( (boost::is_same< void, result_type> :: value) );
          //! Works with boost::result_of.
          template <class Seq>
            struct result
            {
                typedef void type;
            };

        public:
          //! Metafunction for unfused types.
          template< class T_FUNCTOR > struct sequence 
          {
            //! Result type
            typedef typename boost::mpl::push_front
                    <
                      typename boost::function_types::parameter_types
                               < T_FUNCTOR >  :: type,
                      const t_Key&
                    > :: type type; 
          };

          //! Constructor.
          XmlFusedFactory() : node_factory_( new t_NodeFactory ),
                          attribute_factory_( new t_AttributeFactory ),
                          value_factory_( new t_ValueFactory ) {}
          //! Copy Constructor.
          XmlFusedFactory   ( const XmlFusedFactory& _c ) 
                      : node_factory_( _c.node_factory_ ),
                        attribute_factory_( _c.attribute_factory_ ),
                        value_factory_( _c.value_factory_ ) {}

          //! Searches amongst all factories for a function key and creates it.
          template< class T_SEQUENCE > 
            void operator()( T_SEQUENCE& _seq ) const
            {
              BOOST_STATIC_ASSERT(    is_node_sequence<T_SEQUENCE> :: value
                                   or is_attribute_sequence< T_SEQUENCE > :: value 
                                   or is_value_sequence<T_SEQUENCE> :: value );
              operator_( _seq );
            }

          //! Connects a functor to the node factory.
          template< class T_FUNCTOR >
            t_NodeReturn connect_node( const t_Key& _key, const T_FUNCTOR &_functor )
              { return connect( _key, "", _functor ); }
          //! Connects a functor to the node factory.
          template< class T_FUNCTOR >
            t_NodeReturn connect_node( const t_Key& _key, 
                                       const t_Help& _help, 
                                       const T_FUNCTOR &_functor )
            {
              node_factory_->connect( _key, _help, _functor );
              return t_NodeReturn(*this);
            }
          //! Connects a functor to the attribute factory.
          template< class T_FUNCTOR >
            t_AttributeReturn connect_attribute( const t_Key& _key, const T_FUNCTOR &_functor )
              { connect( _key, "", _functor ); }
          //! Connects a functor to the attribute factory.
          template< class T_FUNCTOR >
            t_AttributeReturn connect_attribute( const t_Key& _key,
                                                const t_Help& _help, 
                                                const T_FUNCTOR &_functor )
            { 
              attribute_factory_->connect( _key, _help, _functor );
              return t_AttributeReturn(*this);
            }
          //! Connects a functor to the value factory.
          template< class T_FUNCTOR >
            t_ValueReturn connect_value( const t_Key& _key, const T_FUNCTOR &_functor )
              { return connect( _key, "", _functor ); }
          //! Connects a functor to the value factory.
          template< class T_FUNCTOR >
            t_ValueReturn connect_value( const t_Key& _key,
                                         const t_Help& _help, 
                                         const T_FUNCTOR &_functor )
            {
              value_factory_->connect( _key, _help, _functor );
              return t_ValueReturn(*this); 
            }
               
          //! Returns true if factory key exists.
          bool exists( const t_Key& _key )
            { return     node_factory_->exists( _key ) 
                     or attribute_factory_->exists( _key )
                     or value_factory_->exists( _key ); }

          //! Helper class to check conversion to boost::function<t_Continuator>
          template< class T_FUNCTOR >
            struct is_node : public boost::is_convertible
                      < T_FUNCTOR, boost::function< t_NodeFunction > > {};
          //! Helper class to check conversion to boost::function<t_Updater>
          template< class T_FUNCTOR >
            struct is_attribute : public boost::is_convertible
                      < T_FUNCTOR, boost::function< t_AttributeFunction > > {};
          //! Helper class to check conversion to boost::function<t_Statistics>
          template< class T_FUNCTOR >
            struct is_value : public boost::is_convertible
                      < T_FUNCTOR, boost::function< t_ValueFunction > > {};

        private:
          //! Helper class to check conversion to boost::function<t_Continuator>
          template< class T_SEQUENCE >
            struct is_node_sequence : public boost::mpl::equal 
                      <
                        T_SEQUENCE,
                        typename sequence< t_NodeFunction > :: type,
                        boost::is_same<boost::mpl::_1,boost::mpl::_2>
                      >  {};
          //! Helper class to check conversion to boost::function<t_Updater>
          template< class T_SEQUENCE >
            struct is_attribute_sequence : public boost::mpl::equal
                      <
                        T_SEQUENCE,
                        typename sequence< t_AttributeFunction > :: type,
                        boost::is_same< boost::mpl::_1,  boost::mpl::_2 >
                      > :: type {};
          //! Helper class to check conversion to boost::function<t_Statistics>
          template< class T_SEQUENCE >
            struct is_value_sequence : public boost::mpl::equal
                      <
                        T_SEQUENCE,
                        typename sequence< t_ValueFunction > :: type,
                        boost::is_same< boost::mpl::_1,  boost::mpl::_2 >
                      > :: type {};

          //! Calls a node function.
          template< class T_SEQUENCE >
            typename boost :: enable_if< is_node_sequence<T_SEQUENCE> >::type
              operator_( T_SEQUENCE& _seq ) const { (*node_factory_)( _seq ); }
          //! Calls an attribute function.
          template< class T_SEQUENCE >
            typename boost :: enable_if< is_attribute_sequence<T_SEQUENCE> >::type
              operator_( T_SEQUENCE& _seq ) const { (*attribute_factory_)( _seq ); }
          //! Calls a value function.
          template< class T_SEQUENCE >
            typename boost :: enable_if< is_value_sequence<T_SEQUENCE> >::type
              operator_( T_SEQUENCE& _seq ) const { (*value_factory_)( _seq ); }

          //! Type of the node factory.
          typedef ::LaDa::Factory::FusedFactory< t_NodeFunction, t_Key > t_NodeFactory;
          //! Type of the node factory.
          typedef ::LaDa::Factory::FusedFactory< t_AttributeFunction, t_Key > t_AttributeFactory;
          //! Type of the node factory.
          typedef ::LaDa::Factory::FusedFactory< t_ValueFunction, t_Key > t_ValueFactory;

          //! The node factory.
          boost::shared_ptr<t_NodeFactory> node_factory_;
          //! The attribute factory.
          boost::shared_ptr<t_AttributeFactory> attribute_factory_;
          //! The value factory.
          boost::shared_ptr<t_ValueFactory> value_factory_;
      };

    template< class T_FUNCTION, class T_KEY > 
      std::ostream& operator<<( std::ostream& _stream, 
                                const XmlFusedFactory<T_FUNCTION, T_KEY>& _factory )
      {
        _stream << *_factory.node_factory_ << "\n"
                << *_factory.attribute_factory_ << "\n"
                << *_factory.value_factory_;
        return _stream;
      }

  } // namespace Factory
} // namespace LaDa

#endif 
