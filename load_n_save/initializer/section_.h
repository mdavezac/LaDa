//
//  Version: $Id: section_.h 1231 2009-07-17 05:12:39Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_SECTION_H_
#define _LADA_LOADNSAVE_INITIALIZER_SECTION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/utility/result_of.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/fusion/include/fold.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/static_assert.hpp>

#include "../tree/tree.h"

#include "../grammar/grammar.h"
#include "../dynamic_expression.h"
#include "../transforms/getsections.h"
#include "../transforms/getname.h"
#include "../transforms/gettag.h"
#include "../transforms/getcontent.h"
#include "tagged_options.h"
#include "tagged_sections.h"
#include "option_.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      namespace fusion = boost::fusion;
      class Section_
      {
          //! Type of the xml-like tree.
          typedef tree::Base t_Tree;
     
          //! Returns with error codes and such.
          struct FlowControl;
#         include "control_result_section.impl.h"

        public:

          //! Constructor.
          Section_( t_Tree const& _tree ) : tree_(_tree) {}
          //! Copy Constructor.
          Section_( const Section_& _c ) : tree_(_c.tree_), flow_control(_c.flow_control) {}
         
          //! Result type of the functor.
          typedef bool result_type;
     
          template<class T>
            result_type operator&( T const& _t ) const { return (*this)(_t); }
          //! Do nothing
          result_type operator()( boost::mpl::bool_<false> const& ) const { return true; } 

          //! Functor.
          template< class  T_EXPR > 
            typename boost::disable_if
            <
              typename details::is_view< T_EXPR > :: type,
              result_type
            > :: type operator()( T_EXPR const& _expr ) const;

          //! Goes here when inside an optional group (of sections) with more than one section.
          template< class  T_EXPR >
            typename boost::enable_if
            <
              typename details::is_view< T_EXPR > :: type,
              result_type
            > :: type operator()( T_EXPR const& _expr ) const
              { return opt::while_true( _expr, *this ); }
           //! Goes here to iterate over different optional groups of sections.
           template< class  T_EXPR >
             bool operator()( transform::OptionalGroups<T_EXPR const> const& _expr ) const;
          
          //! Calls virtual function of dynamic expression \a _expr.
          result_type operator()( Expression const &_expr ) const { return _expr(*this); }
#         ifdef DSEL_TYPENAME
#           error DSEL_TYPENAME already exists.
#         endif
#         define DSEL_TYPENAME                                         \
               boost::proto::expr                                      \
               <                                                       \
                 boost::proto::tag::terminal,                          \
                 boost::proto::term< grammar::details::external_type<T_TYPE> > \
               > 
          //! initializes using a dynamic expression \a _expr.
          template< class T_TYPE >
            bool operator()( grammar::Dsel<DSEL_TYPENAME> const &_expr ) const
              { return _expr(*this); }
#         undef DSEL_TYPENAME
     
        protected:
          //! Prints unknown options. 
          void print_unknowns_( std::string const& _name, tree::Section const & _sec ) const;
          //! Does actual parsing. No grammar check except for option values.
          template< class T_EXPR >
            void parse_( T_EXPR const &_expr,
                         tree::Section const& _section,
                         FlowControl & _c ) const;
          //! Does nothing
          void parse_( boost::mpl::bool_<false> const &,
                       tree::Section const&,
                       FlowControl &_c ) const;
          //! The current text.
          t_Tree const& tree_;
          //! Holds a FlowControl instance.
          FlowControl flow_control;
     
      };

#     include "control_result_section.impl.h"

      template< class  T_EXPR >
        bool Section_ :: operator()( transform::OptionalGroups<T_EXPR const> const& _expr ) const
        {
          namespace bl = boost::lambda;
          Section_ other(*this);
          other.flow_control.quiet = true;
//         other.flow_control.check_unknown_sections = false;
//         other.flow_control.check_unknown_options = false;
          return not opt::while_true( _expr.view, not bl::bind( other, bl::_1 ) );
        }

      template< class  T_EXPR >
        typename boost::disable_if
        <
          typename details::is_view< T_EXPR > :: type,
          Section_::result_type
        > :: type Section_::operator()( T_EXPR const& _expr ) const
        {
          BOOST_STATIC_ASSERT
            (( boost::proto::matches<T_EXPR, grammar::OrnateSection>::type::value ));
        
          FlowControl result(flow_control);
          std::string const name( transform::GetName()(_expr) );
          size_t tags( transform::GetTag()( _expr ) );
          result.is_required = tags & tags::section::required;
          result.is_in_xml_tree = false;

          t_Tree :: const_iterator :: subsection i_first( tree_.subsections_begin(name) );
          t_Tree :: const_iterator :: subsection const i_last( tree_.subsections_end(name) );
          if( i_first == i_last ) return result.result(_expr);

          typedef typename boost::result_of< transform::GetContent(T_EXPR) >::type t_Content;
          t_Content const content( transform::GetContent()(_expr) );
          TaggedOptions req_options( content, tags::option::required );
          TaggedSections req_sections( content, tags::section::required );
          for(; i_first != i_last; ++i_first )
          {
            if( not i_first->fresh ) continue;
            if( not id_options( *i_first, content ) ) continue;

            result.is_in_xml_tree = true;
        
            // Found a section which matches.
            // Now checks it is valid. 
            if( not result.do_parse(_expr, *i_first) ) continue;
            i_first->fresh = false; // marks section as read.
        
            // finally parses.
            parse_( transform::GetContent()(_expr), *i_first, result );
            result.did_parse( name );
            print_unknowns_( name, *i_first );
            break;
          }
          return result.result( _expr );
        }

      
      template< class T_EXPR >
        void Section_::parse_( T_EXPR const &_expr,
                               tree::Section const& _section,
                               FlowControl &_result ) const
        {
          BOOST_STATIC_ASSERT
            (( boost::proto::matches<T_EXPR, grammar::Content>::type::value ));
          namespace bl = boost::lambda;

          // parse options.
          _result.parsed_options = true;
          fusion::for_each
          ( 
            transform::GetOptions()(_expr),
            bl::var(_result.parsed_options) &= bl::bind(Option_(_section, false), bl::_1)
          );
     
          // parse sections.
          _result.parsed_subsections = ( Section_(_section) )(transform::GetSections()(_expr));
        }

      inline void Section_::parse_( boost::mpl::bool_<false> const &,
                                    tree::Section const&,
                                    FlowControl &_c ) const 
      {
        _c.parsed_options = true;
        _c.parsed_subsections = true;
      }

    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa


#endif
