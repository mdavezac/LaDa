//
//  Version: $Id: option_.h 1227 2009-07-14 02:17:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_OPTION_H_
#define _LADA_LOADNSAVE_INITIALIZER_OPTION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/utility/result_of.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/fusion/include/is_view.hpp>
#include <boost/static_assert.hpp>

#include "opt/while_true.h"

#include "../tree/tree.h"

#include "../grammar/grammar.h"
#include "../transforms/getname.h"
#include "../transforms/gettag.h"
#include "../transforms/getaction.h"
#include "../transforms/getvalues.h"
#include "../transforms/getdefault.h"
#include "../transforms/group_optionals.h"
#include "value_.h"

#include <opt/debug.h>

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      namespace details
      {
        //! General is_view metafunctor.
        template< class T_EXPR >
          class is_view : public boost::fusion::traits::is_view<T_EXPR> {};
        //! Specialization when instanciated with boost::mpl::false_.
        template<>
          class is_view<boost::mpl::false_> : public boost::mpl::false_ {};
      };

      namespace fusion = boost::fusion;
      class Option_
      {
          //! Type of the xml-like tree.
          typedef tree::Section t_Tree;
          //! A structure to control the parsing flow, results, and errors.
          struct FlowControl
          {
            //! Controls whether to match or parse-and-match.
            bool match_only;
            //! Controls which tags to match.
            size_t match_tag;
            //! Whether can have many of this option.
            bool is_many;
            //! Whether is an id option.
            bool is_id;
            //! Whether is an id option.
            bool is_required;
            //! Whether this option has been found.
            bool found;
            //! Whether this option has been found more than once.
            bool found_many;
            //! Whether the option exists in the xml tree.
            bool is_in_xml_tree;
            //! Whether to perform defaults when no match is found.
            bool no_default;
          
            //! Default Constructor.
            FlowControl();
            //! Copy Constructor.
            FlowControl   ( FlowControl const& _c );
            //! True if matches correct tag
            bool good_tag( size_t _tag ) const;
            //! When a new matching option is found.
            void found_another() { found ? found_many = true: found = true; }
            //! Re-inits results.
            void init( size_t _tag );
            //! Return true/false, + prints error messages.
            bool operator()( std::string const& ) const;
            //! Whether to check for default value.
            bool do_default() const { return not ( no_default or is_id or is_required or found ); }

          };

        public:
          //! Constructor.
          Option_   ( t_Tree const& _tree, bool _is_matching,
                      size_t _tag = tags::option::default_ )
                  : tree_(_tree)
          {
            flow_control.match_only = _is_matching;
            flow_control.match_tag = _tag;
          }
          //! Copy Constructor.
          Option_    ( const Option_& _c )
                   : tree_(_c.tree_), flow_control( _c.flow_control ) {}
         
          //! Result type of the functor.
          typedef bool result_type;
     
          //! Functor.
          template< class  T_EXPR > 
            typename boost::disable_if
            <
              typename details::is_view< T_EXPR > :: type,
              bool
            > :: type operator()( T_EXPR const& _expr ) const;
          //! Goes here when inside an optional group (of options) with more than one option.
          template< class  T_EXPR >
            typename boost::enable_if
            <
              typename details::is_view< T_EXPR > :: type,
              bool
            > :: type operator()( T_EXPR const& _expr ) const
              { return opt::while_true( _expr, *this ); }
           //! Goes here to iterate over different optional groups of options.
           template< class  T_EXPR >
             bool operator()( transform::OptionalGroups<T_EXPR const> const& _expr ) const;
          
          //! Do nothing
          bool operator()( boost::mpl::bool_<false> const& ) const { return true; } 
     
        protected:
          //! The current text.
          t_Tree const& tree_;
          //! Controls parsing flow, results, and errors.
          FlowControl flow_control;
      };

      template< class  T_EXPR >
        bool Option_ :: operator()( transform::OptionalGroups<T_EXPR const> const& _expr ) const
        {
          namespace bl = boost::lambda;
          Option_ other(*this);
          other.flow_control.no_default = true;
          other.flow_control.match_only = true;
          other.flow_control.found = not opt::while_true( _expr.view, not bl::bind( other, bl::_1 ) );
          if( other.flow_control.found or flow_control.match_only ) return true;
          if( flow_control.is_id ) return false;
          if( flow_control.is_required ) return false;
          other.flow_control.no_default = false;
          return not opt::while_true( _expr.view, not bl::bind( other, bl::_1 ) );
        }

      template< class  T_EXPR >
        typename boost::disable_if
        <
          typename details::is_view< T_EXPR > :: type,
          bool
        > :: type Option_ :: operator()( T_EXPR const& _expr ) const
        {
          BOOST_STATIC_ASSERT
            (( boost::proto::matches<T_EXPR, grammar::OrnateOption>::type::value ));
          
          std::string const name( transform::GetName()(_expr) );
          std::string const section_name( tree_.name );
          size_t const tags( transform::GetTag()(_expr) );
          FlowControl result( flow_control );

          result.init( tags );
          if( not result.good_tag( tags ) ) return true;
        
 
          t_Tree :: const_iterator :: option i_first( tree_.options_begin(name) );
          t_Tree :: const_iterator :: option const i_last( tree_.options_end(name) );
 
          bool found(false);
          typedef typename boost::result_of< transform::GetAction(T_EXPR) > :: type t_Action;
          t_Action const action( transform::GetAction()(_expr) );
          Value_ const value_(section_name);
          result.is_in_xml_tree = i_first != i_last;
          for(; i_first != i_last; ++i_first )
          {
            if( not i_first->fresh ) continue;
            if
            ( 
             not value_( transform::GetValues()(_expr), action,
                         i_first->name, i_first->value) 
            ) continue;

            result.found_another(); 
            i_first->fresh = false;
          }

          const bool has_default
          (
            not boost::is_same
              < 
                typename boost::result_of< transform::GetDefault(T_EXPR) > :: type, 
                boost::mpl::bool_<false> 
              > :: type :: value
          );
          LADA_ASSERT( has_default and tags & tags::option::required, "This is a bug.\n");
          LADA_ASSERT( has_default and tags & tags::option::id, "This is a bug.\n");

          if( result.do_default() and  has_default )
          {
            std::string const default_( convert_value_to_regex( transform::GetDefault()(_expr) ) );
            const bool good
            (
              StraightValue_()( default_, boost::mpl::bool_<false>(), action )
            );
            LADA_ASSERT( not good,   "Bug: could not parse default value of option " 
                                   + name + " in section " + section_name + ".\n" );
            result.found = true;
          }
          return result( name );
        }
      
    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa


#endif
