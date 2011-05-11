//
//  Version: $Id: find_named_argument.h 1211 2009-07-04 02:03:24Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_FIND_NAMED_ARGUMENT_H_
#define _LADA_LOADNSAVE_TRANSFORMS_FIND_NAMED_ARGUMENT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/mpl/less_equal.hpp>
#include <boost/mpl/less.hpp>

#include "../grammar/option_terminal.h"
#include "../grammar/section_terminal.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;

      namespace details
      {

        template< size_t N, class Name, class Default >
          struct ReturnValue : public proto::if_
            <
              boost::mpl::less< boost::mpl::long_<N>, proto::arity_of<proto::_> >(),
              proto::if_
              < 
                proto::matches< proto::_child_c<N>, Name >(), 
                proto::_value( proto::_right(proto::_child_c<N>) ),
                ReturnValue<N+1,Name, Default>()
              >(),
              Default()
            >  {};
      }

      template< class Name, class Default>
        struct GetNamedArg : public proto::or_
          <
            proto::when
            < 
              grammar::option_function, 
              details::ReturnValue< 2, Name, Default >() 
            >,
            proto::when
            < 
              grammar::section_function, 
              details::ReturnValue< 2, Name, Default >() 
            >,
            proto::when
            <
              proto::assign< grammar::section_function, grammar::Content >,
              GetNamedArg<Name, Default>( proto::_left )
            >,
            proto::when
            <
              proto::assign< grammar::option_function, grammar::Values >,
              GetNamedArg<Name, Default>( proto::_left )
            >,
            proto::otherwise< Default() >
          > {};

    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif
