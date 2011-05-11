//
//  Version: $Id: dynamic_expression.h 1201 2009-06-22 05:26:09Z davezac $
//

#ifndef _LADA_LOADNSAVE_DYNAMIC_EXPRESSION_H_
#define _LADA_LOADNSAVE_DYNAMIC_EXPRESSION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/shared_ptr.hpp>
#include <boost/proto/proto.hpp>

#include "expression_tags.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace initializer
    {
      class Initializer;
      class Section_;
    };

    namespace details
    {
      //! Abstracts out a type  via derived details::TypedExpression.
      struct ExpressionVirtualBase
      {
        //! Destructor.
        virtual ~ExpressionVirtualBase() {};
        //! call back.
        virtual bool operator()( initializer::Initializer const& _a ) const = 0;
        //! call back.
        virtual bool operator()( initializer::Section_ const& _a ) const = 0;
      };


      //! Type of expression is abstracted out and returned via a callback.
      template< class T_EXPR >
        struct TypedExpression : public ExpressionVirtualBase
        {
          //! Type of the expression.
          typedef T_EXPR t_Expression;
      
          //! Constructor.
          TypedExpression( t_Expression const &_expr ) : expr_(_expr) {};
          //! Copy Constructor.
          TypedExpression( TypedExpression const &_c ) : expr_(_c.expr_) {};
          //! Destructor.
          virtual ~TypedExpression() {}
          //! call back.
          virtual bool operator()( initializer::Initializer const& _a ) const { return _a & expr_; }
          //! call back.
          virtual bool operator()( initializer::Section_ const& _a ) const
            { return call_(_a, typename expression_types::is_section<t_Expression>::type() ); }

          private:
            //! Tag dispatched call to something.
            template<class T_CALL>
              bool call_( T_CALL const&_a, boost::mpl::true_ const& ) const
                { return _a(expr_); }
            //! Tag dispatched non-call to something.
            template<class T_NOCALL>
              bool call_( T_NOCALL const&_a, boost::mpl::false_ const& ) const
                { return true; }
            //! Holds the abstracted type.
            t_Expression const& expr_;
        };

    }

    namespace grammar
    {
      //! Specialization for dynamic expression.
      template<>  struct Dsel<dynamic_expression::type>
        : public proto::extends< dynamic_expression::type, Dsel<dynamic_expression::type>, Domain >
        {
          protected:
            //! The base type.
            typedef proto::extends< dynamic_expression::type, Dsel<dynamic_expression::type>, Domain > t_Base;
            //! Holds a pointer to the implementation.
            boost::shared_ptr< load_n_save::details::ExpressionVirtualBase > impl_;

          public:
            //! A general constructor to abstract out expression types.
            template< class T_EXPR > 
              Dsel( T_EXPR const& _expr ) : impl_( new load_n_save::details::TypedExpression<T_EXPR>(_expr) ) {}
            //! Copy constructor.
            Dsel( Dsel const& _a ) : t_Base(_a), impl_(_a.impl_) {};
            //! call back.
            bool operator()( initializer::Initializer const& _a ) const { return (*impl_)(_a); }
            //! call back.
            bool operator()( initializer::Section_ const& _a ) const { return (*impl_)(_a); }

            //! Resets to another static expression.
            template< class T_EXPR >
              void operator=( T_EXPR const &_a ) const
              { 
                boost::shared_ptr< load_n_save::details::ExpressionVirtualBase > dummy(_a);
                impl_.swap( dummy ); 
              }
            //! Copies (shallow) another dynamic expression.
            void operator=( Dsel<dynamic_expression::type> const &_a ) 
            {
              boost::shared_ptr< load_n_save::details::ExpressionVirtualBase > dummy(_a.impl_);
              impl_.swap( dummy ); 
            }
        };
    }

    //! Typedef to a dynamic expression.
    typedef grammar::Dsel<grammar::dynamic_expression::type> Expression;
 
  } // namespace load_n_save

} // namespace LaDa


#endif
