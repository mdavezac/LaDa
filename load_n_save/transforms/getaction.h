//
//  Version: $Id: getaction.h 1212 2009-07-04 04:36:19Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_GETACTION_H_
#define _LADA_LOADNSAVE_TRANSFORMS_GETACTION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "find_named_argument.h"
#include "false.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;

      namespace details
      {

        template< size_t N >
          struct ReturnAction : public proto::if_
            <
              boost::mpl::less< boost::mpl::long_<N>, proto::arity_of<proto::_> >(),
              proto::if_
              < 
                proto::matches< proto::_child_c<N>, grammar::action_grammar >(), 
                var( proto::_right(proto::_child_c<N>) ),
                ReturnAction<N+1>()
              >(),
              boost::mpl::bool_<false>()
            >  {};
      }
      //! Returns action of section or option.
      struct GetAction  : public proto::or_
          <
            proto::when
            < 
              grammar::option_function, 
              details::ReturnAction<2>()
            >,
            proto::when
            < 
              grammar::section_function, 
              details::ReturnAction<2>()
            >,
            proto::when
            <
              proto::assign< grammar::option_function, grammar::TypedValues<proto::_> >,
              details::ReturnAction<2>( proto::_left )
            >,
            proto::otherwise< boost::mpl::bool_<false>() >
          > {};
  
    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif
