//
//  Version: $Id: false.h 1122 2009-05-18 02:27:54Z Mayeul $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_FALSE_H_
#define _LADA_LOADNSAVE_TRANSFORMS_FALSE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/tags.hpp>

#include <boost/fusion/include/single_view.hpp>

#include <boost/mpl/bool.hpp>

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      //! Placeholder when subsection, action, or what not, is not present.
      struct false_ : public boost::proto::when
       <
         boost::proto::_,
         boost::fusion::single_view< boost::mpl::bool_<false> >()
       > {};


    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif
