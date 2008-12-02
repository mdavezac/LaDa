//
//  Version: $Id: factory.h 860 2008-11-17 18:37:10Z davezac $
//
#ifndef _LADA_FACTORY_XMLFACTORY_H_
#define _LADA_FACTORY_XMLFACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <map>

#include <boost/fusion/sequence/intrinsic/at.hpp>


#include "xmlfusedfactory.h"

namespace LaDa
{
  //! Holds factory related objects.
  namespace Factory
  {
    //! \brief Intermediate instanciation to make sure all base classes in
    //!        XmlFactory are correctly initialized.
    template
    < 
      class T_FUNCTION, 
      class T_KEY = XmlKey,
      class SIZE = typename boost::fusion::result_of::size
                   <
                     typename boost::function_types::parameter_types< T_FUNCTION > :: type
                   > :: type
    > class XmlFactory;

    //! An Xml factory for function taking no argument.
    template< class T_FUNCTION, class T_KEY>
      class XmlFactory<T_FUNCTION, T_KEY, boost::mpl::long_<0> >  
        : public XmlFusedFactory< T_FUNCTION, T_KEY >
      {
          //! Base class.
          typedef XmlFusedFactory<T_FUNCTION, T_KEY> t_Base;
          //! Type of the return.
          typedef typename t_Base :: result_type  result_type;
          //! Type of the Xmlnode.
          typedef typename t_Base :: t_XmlElement  t_XmlElement;
          //! Type of the Xml attribute.
          typedef typename t_Base :: t_XmlAttribute  t_XmlAttribute;
          //! Parameter list
          typedef typename boost::function_types::parameter_types< T_FUNCTION > :: type t_Parameters;

        public:
          //! Type of the key.
          typedef typename t_Base :: t_Key  t_Key;
          //! Type of the key.
          typedef typename t_Base :: t_Function  t_Function;

          result_type operator()(const t_Key& _key, const t_XmlElement& _node ) const
          {
            boost::fusion::vector< const t_Key&, const t_XmlElement& > vector( _key, _node );
            return t_Base::operator()( vector );
          }
          result_type operator()(const t_Key& _key, const t_XmlAttribute& _attribute ) const
          {
            boost::fusion::vector< const t_Key&, const std::string& > vector( _key, _attribute );
            return t_Base::operator()( vector );
          }
          result_type operator()(const t_Key& _key) const
          {
            boost::fusion::vector< const t_Key& > vector( _key );
            return t_Base::operator()( vector );
          }

          using t_Base::connect_node;
          using t_Base::connect_attribute;
          using t_Base::connect_value;
      };

    //! An Xml factory for function taking a single argument.
    template< class T_FUNCTION, class T_KEY>
      class XmlFactory<T_FUNCTION, T_KEY, boost::mpl::long_<1> >  
        : public XmlFusedFactory< T_FUNCTION, T_KEY >
      {
          //! Base class.
          typedef XmlFusedFactory<T_FUNCTION, T_KEY> t_Base;
          //! Type of the return.
          typedef typename t_Base :: result_type  result_type;
          //! Type of the Xmlnode.
          typedef typename t_Base :: t_XmlElement  t_XmlElement;
          //! Type of the Xml attribute.
          typedef typename t_Base :: t_XmlAttribute  t_XmlAttribute;
          //! Parameter list
          typedef typename boost::function_types::parameter_types< T_FUNCTION > :: type t_Parameters;
          //! First Argument
          typedef typename boost::fusion::result_of::at
                           < 
                             t_Parameters, 
                             boost::mpl::int_<0> 
                           > :: type t_Arg1;

        public:
          //! Type of the key.
          typedef typename t_Base :: t_Key  t_Key;
          //! Type of the key.
          typedef typename t_Base :: t_Function  t_Function;

          result_type operator()(const t_Key& _key, t_Arg1 _arg1, const t_XmlElement& _node ) const
          {
            boost::fusion::vector< const t_Key&, t_Arg1, const t_XmlElement& > vector( _key, _arg1, _node );
            return t_Base::operator()( vector );
          }
          result_type operator()(const t_Key& _key, t_Arg1 _arg1, const t_XmlAttribute& _attribute ) const
          {
            boost::fusion::vector< const t_Key&, t_Arg1, const t_XmlAttribute& >
              vector( _key, _arg1, _attribute );
            return t_Base::operator()( vector );
          }
          result_type operator()(const t_Key& _key, t_Arg1 _arg1 ) const
          {
            boost::fusion::vector< const t_Key&, t_Arg1 > vector( _key, _arg1 );
            return t_Base::operator()( vector );
          }
          using t_Base::connect_node;
          using t_Base::connect_attribute;
          using t_Base::connect_value;
      };

  }
}

#endif 
