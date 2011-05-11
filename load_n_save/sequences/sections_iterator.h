//
//  Version: $Id: sections_iterator.h 1196 2009-06-20 17:14:26Z davezac $
//

#ifndef _LADA_LOADNSAVE_STATIC_EXPRESSION_H_
#define _LADA_LOADNSAVE_STATIC_EXPRESSION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/mpl/bool.hpp>
#include <boost/fusion/include/iterator_facade.hpp>

#include "../expression_tags.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace sequences
    {
      template<class T_EXPR>
        struct SectionsIterator
          : public boost::fusion::iterator_facade
                   <
                     SectionsIterator<T_EXPR>,
                     boost::fusion::forward_traversal_tag
                   >
        {
          //! Type of the base.
          typedef boost::fusion::iterator_facade
                  <
                    SectionsIterator<T_EXPR>,
                    boost::fusion::forward_traversal_tag
                  > t_Base;
          //! Returns first element and rest of expression.
          struct first_and_others;
          public:
            //! Type of the expression.
            typedef T_EXPR t_Expression;

            //! Constructor.
            Sections( t_Expression const& _expr )
        };

      template<class T_EXPR>
        struct SectionsIterator::first_and_others
          : public proto::or_
          <
            proto::when
            < 
              grammar::OrnateSection, 
              join( protect( proto::_ ), false_() )
            >,
            proto::when
            <
              proto::or_
              <
                proto::plus< grammar::OrnateOptions, grammar::OrnateSection >,
                proto::plus< grammar::Content, grammar::OrnateSection >,
                proto::plus< grammar::OrnateSections, grammar::OrnateSection >,
              >,
              join( GetSections(proto::_right), protect(proto::_left) )
            >
            proto::otherwise< join( false_(), false_() ) >
          > {};

    } // namespace sequences.
  } // namespace load_n_save
} // namespace LaDa

