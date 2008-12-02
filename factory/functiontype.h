//
//  Version: $Id: factory.h 860 2008-11-17 18:37:10Z davezac $
//
#ifndef _LADA_FACTORY_FUNCTIONTYPE_H_
#define _LADA_FACTORY_FUNCTIONTYPE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/identity.hpp>

#include <boost/fusion/adapted/mpl.hpp>

#include <boost/function_types/result_type.hpp>
#include <boost/function_types/components.hpp>
#include <boost/function_types/function_type.hpp>

namespace LaDa
{
  //! Holds factory related objects.
  namespace Factory
  {

    //! Adds a parameter to function parameter list.
    template< class T_FUNCTION, class T_NEWARG >
      struct push_back_arg : public boost::function_types::function_type
                             < 
                               typename boost::mpl::push_back
                               < 
                                 typename boost::function_types::components
                                 <
                                   T_FUNCTION, 
                                   boost::mpl::identity<boost::mpl::_1>
                                 >::type,
                                 T_NEWARG
                               > :: type
                             >  {};

  } // namespace Factory.
} // LaDa



#endif

