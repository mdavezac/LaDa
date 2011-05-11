//
//  Version: $Id: sections.h 1231 2009-07-17 05:12:39Z davezac $
//

#ifndef _LADA_LOADNSAVE_SECTION_DEFS_H_
#define _LADA_LOADNSAVE_SECTION_DEFS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _LADADEBUG
# include <iostream>
#endif
#include <boost/fusion/include/at_c.hpp>

#include <boost/proto/expr.hpp>
#include <boost/proto/extends.hpp>
#include <boost/proto/tags.hpp>
#include <boost/proto/traits.hpp>
#include <boost/proto/operators.hpp>
#include <boost/proto/domain.hpp>

#include <opt/types.h>

#include "options.h"
#include "section_terminal.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace grammar
    {
      namespace proto = boost::proto;
     
      // forward declaration.
      struct OrnateSections;
      
   
      //! \details Assignement for an option. Partial.
      //! \brief This grammar works when the assignee of a section include options.
      //!        The case when they don't is dealt with in OrnateSection.
      //!        This way, options must come at the beginning.
       struct Content : public proto::or_
         <
           OrnateOptions,
           OrnateSections,
           proto::plus< OrnateOptions, OrnateSections >,
           proto::plus< Content, OrnateSections >
         > {};

      //! Sections can ornated with (assigned) options and (sub-)sections using plus.
      struct OrnateSection : public proto::or_
        <
          section_function,
          dynamic_expression,
          external_type,
          action_type,
          proto::assign< section_function, Content >
        > {};
   
      //! Sections of equivalent depth are linked via plus.
      struct OrnateSections : public proto::or_
        <
          OrnateSection,
          proto::plus< OrnateSections, OrnateSections >,
          proto::logical_or< OrnateSections, OrnateSections >
        > {};
    } // namespace grammar

  } // namespace load_n_save

} // namespace LaDa


#endif
