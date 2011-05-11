//
//  Version: $Id: expression_tags.h 1198 2009-06-21 21:37:05Z davezac $
//

#ifndef _LADA_LOADNSAVE_EXPRESSION_TYPES_H_
#define _LADA_LOADNSAVE_EXPRESSION_TYPES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "grammar/sections.h" 
#include "grammar/options.h" 

namespace LaDa 
{
  namespace load_n_save
  {
    namespace expression_types
    {
      //! An option tag.
      struct option_tag {};
      //! A section tag.
      struct section_tag {};
      //! An "other" tag.
      struct other_tag {};


      //! True is T_EXPR is a section.
      template< class T_EXPR >
        struct is_section : public boost::proto::matches 
              < 
                T_EXPR, 
                grammar::OrnateSection
              > :: type {};
      //! True is T_EXPR is an option.
      template< class T_EXPR >
        struct is_option : public boost::proto::matches 
              < 
                T_EXPR, 
                grammar::OrnateOption
              > :: type {};
      //! True is T_EXPR is an "other".
      template< class T_EXPR >
        struct is_other : public boost::mpl::and_
              < 
                typename boost::mpl::not_< is_section<T_EXPR> > :: type,
                typename boost::mpl::not_< is_option<T_EXPR> > :: type
              > :: type {};
      template< class T_EXPR, class ENIF = void >
        struct tag
        { 
          //! default tag is other.
          typedef other_tag type; 
        };

      //! Matches a section.
      template< class T_EXPR >
        struct tag
        < 
          T_EXPR,
          typename boost::enable_if<typename is_section<T_EXPR>::type>::type
        > 
        {
          //! Type is a section.
          typedef section_tag type;
        };

      //! Matches an option.
      template< class T_EXPR >
        struct tag
        < 
          T_EXPR,
          typename boost::enable_if<typename is_option<T_EXPR>::type>::type
        >  
        {
          //! Type is a section.
          typedef option_tag type;
        };

    }
  } // namespace load_n_save

} // namespace LaDa


#endif
