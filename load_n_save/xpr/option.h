//
//  Version: $Id: option.h 1268 2009-08-12 04:56:27Z davezac $
//

#ifndef _LADA_LNS_XPR_TREE_OPTION_H_
#define _LADA_LNS_XPR_TREE_OPTION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/utility/enable_if.hpp>

#include "../string_type.h"
#include "../action/action.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      //! An option/value pair.
      class Option
      {
        public:
          //! Name of the option
          t_String name;
          //! Tag.
          size_t tag;
          //! Help string.
          t_String help;

          //! Constructor.
          Option   ( const t_String _name = "" )
                 : name(_name), tag(0), help(""), action_() {}
          //! Copy Constructor.
          Option   ( Option const& _c )
                 : name(_c.name), tag(_c.tag), help(_c.help), action_(_c.action_) {}
          //! calls the action.
          bool action(t_String const& _str) const
            { return action_ ? (*action_)(_str): _str.size() == 0; }
          //! calls the action.
          t_String action() const { return action_ ? (*action_)(): ""; }
          //! Assigns a default value if it exists.
          bool assign_default() const { return action_ ? action_->assign_default(): false; }
          //! Sets the action.
          template< class T0, class T1 >
            typename boost::disable_if< action_::is_special_action<T0> > :: type 
              set_action( T0& _ref, T1 const& _def) 
                { action_.reset( new action_::Action<T0,T1>(_ref, _def) ); }
          //! Sets the action.
          template< class T0, class T1 >
            typename boost::enable_if< action_::is_special_action<T0> > :: type 
              set_action( T0& _ref, T1 const& _def) 
                { action_.reset( new action_::SpecialAction<T0,T1>(_ref, _def) ); }
          //! Deletes action.
          void set_action( boost::mpl::false_, boost::mpl::false_ )
            { boost::shared_ptr<action_::ActionBase>().swap( action_ ); }
        protected:
          //! Action.
          boost::shared_ptr<action_::ActionBase> action_;
      };

    } // namespace parser

  } // namespace load_n_save
} // namespace LaDa

#endif
