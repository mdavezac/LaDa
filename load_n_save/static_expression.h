//
//  Version: $Id: static_expression.h 1195 2009-06-20 00:28:59Z davezac $
//

#ifndef _LADA_LOADNSAVE_STATIC_EXPRESSION_H_
#define _LADA_LOADNSAVE_STATIC_EXPRESSION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "expression_tags.h"
#include "transforms/getname.h"
#include "transforms/getaction.h"
#include "transforms/gettag.h"

namespace LaDa 
{
  namespace load_n_save
  {
    template< class T_EXPR >
      class StaticExpression
      {
        public:
          //! Type of the expression.
          typedef T_EXPR t_Expression;
          //! Type of this expression.
          typedef typedef expression_type :: tag< T_EXPR > :: type tag;

          struct sequences
          {
            //! Type of the sections sequence.
            typedef load_n_save::sequences::Section<T_EXPR> sections;
            //! Type of the option iterators.
            typedef load_n_save::sequences::Options<T_EXPR> options;
          };

          //! Constructor.
          StaticExpression( t_Expression const &_expr ) : expr_(_expr) {};
          //! Copy Constructor.
          StaticExpression( StaticExpression const &_c ) : expr_(_c.expr_) {};

          //! Returns the name, if any, of this expression.
          typename boost::result_of< transform::GetName( const t_Expression& ) > :: type
            name() const { return transform::GetName()( expr_ ); }

          //! Returns the tags, if any, of this expression.
          typename boost::result_of< transform::GetTag( const t_Expression& ) > :: type
            tag() const { return transform::GetTag()( expr_ ); }

          //! Returns the action, if any, of this expression.
          typename boost::result_of< transform::GetAction( const t_Expression& ) > :: type
            action() const { return transform::GetAction()( expr_ ); }

          //! Returns a fusion sequence of sections. 
          sequences::sections sections() const { return sequences::sections( expr_ ); }
          //! Returns a fusion sequence of options. 
          sequences::options options() const { return sequences::options( expr_ ); }

        protected:
          t_Expression const &expr_;
      };
 
  } // namespace load_n_save

} // namespace LaDa


#endif
