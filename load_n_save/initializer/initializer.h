//
//  Version: $Id: initializer.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_INITIALIZER_H_
#define _LADA_LOADNSAVE_INITIALIZER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/static_assert.hpp>

#include "../tree/tree.h"

#include "../grammar/grammar.h"
#include "../transforms/getsections.h"
#include "../dynamic_expression.h"

#include "section_.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      namespace fusion = boost::fusion;
     
      class Initializer 
      {
        public:
          //! Type of the xml-like tree.
          typedef tree::Base t_Tree;
     
          //! Constructor.
          Initializer() {}
          //! Copy constructor.
          Initializer( const Initializer &_c ) : tree_( _c.tree_ ) {}
          
          //! initializes from an xml tree using a static expression \a _expr.
          template< class T_EXPR > bool operator()( t_Tree const &_tree, T_EXPR const &_expr )
           { tree_ = _tree; return (*this) & _expr; } 
          //! initializes using a dynamic expression \a _expr.
          bool operator&( Expression const &_expr ) const
            { return _expr(*this); }
          //! initializes using a static expression \a _expr.
          template< class T_EXPR >
            typename boost::disable_if
            <
              typename boost::proto::matches< T_EXPR, grammar::dynamic_expression > :: type,
              bool
            > :: type operator&( T_EXPR const &_expr ) const;
          template< class T_EXPR >
            typename boost::enable_if
            <
              typename boost::proto::matches< T_EXPR, grammar::dynamic_expression > :: type,
              bool
            > :: type operator&( T_EXPR const &_expr ) const { return _expr(*this); }

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
            bool operator&( grammar::Dsel<DSEL_TYPENAME> const &_expr ) const
              { return _expr(*this); }
#         undef DSEL_TYPENAME
     
          bool is_loading() const { return true; } 
     
        protected:
          //! Whether we are currently parsing.
          bool is_parsing_;
          //! Xml-like tree.
          t_Tree tree_;
      };
     
     
      template< class T_EXPR >
        typename boost::disable_if
        <
          typename boost::proto::matches< T_EXPR, grammar::dynamic_expression > :: type,
          bool
        > :: type Initializer::operator&( T_EXPR const &_expr ) const
        {
          namespace bl = boost::lambda;
          BOOST_STATIC_ASSERT
            (( boost::proto::matches<T_EXPR, grammar::OrnateSections>::type::value ));
          bool result(true);
          fusion::for_each
          ( 
            transform::GetSections()(_expr),
            bl::var(result) &= bl::bind( Section_(tree_), bl::_1)
          );
          return result;
        }
    } // namespace initializer.
  } // namespace load_n_save

} // namespace LaDa


#endif
