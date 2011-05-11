//
//  Version: $Id: static_expression.h 1194 2009-06-19 00:11:56Z davezac $
//

#ifndef _LADA_LOADNSAVE_STATIC_EXPRESSION_H_
#define _LADA_LOADNSAVE_STATIC_EXPRESSION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/mpl/bool.hpp>
#include <boost/fusion/include/sequence_facade.hpp>

#include "../expression_tags.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace sequences
    {
      template<class T_EXPR>
        struct Sections : public boost::fusion::sequence_facade
                          <
                            Sections<T_EXPR>,
                            boost::fusion::forward_traversal_tag,
                            boost::mpl::true_
                          >
        {
          typedef boost::fusion::sequence_facade
                  <
                    Sections<T_EXPR>,
                    boost::fusion::forward_traversal_tag,
                    boost::mpl::true_
                  > t_Base;
          public:
            //! Type of the expression.
            typedef T_EXPR t_Expression;

            //! Constructor.
            Sections( t_Expression const& _expr ) : expr_(_expr) {};
            //! Constructor.
            Sections( Sections const& _c ) : expr_(_c.expr_) {};

            //! Returns the size of the sequence.
            template< class T_SQ > struct size;
            //! Returns the first iterator.
            template< class T_SQ > struct begin;
            //! Returns the second iterator.
            template< class T_SQ > struct end;
          protected:


            //! The expression for which this is a sequence.
            t_Expression const &expr_;
        };

      template<class T_EXPR> template<class T_SQ>
        struct Sections<T_EXPR> :: size
        {
          protected:
            //! Counts the number of occurences.
            struct count : public proto::or_
              <
                proto::when< grammar::OrnateSection, boost::mpl::int_<1>() >,
                proto::when
                <
                  proto::plus< grammar::OrnateOptions, grammar::OrnateSection >,
                  boost::mpl::int_<1>()
                >,
                proto::when
                <
                  proto::plus< grammar::Content, grammar::OrnateSection >,
                  boost::mpl::plus< count(proto::_left), boost::mpl::int_<1> >()
                >,
                proto::when
                <
                  proto::plus< grammar::OrnateSections, grammar::OrnateSection >,
                  boost::mpl::plus< count(proto::_left), boost::mpl::int_<1> >()
                >,
                proto::when
                <
                  proto::plus< grammar::OrnateSections, grammar::OrnateSections >,
                  boost::mpl::plus<count(proto::_left), count(proto::_right)>()
                >,
                proto::otherwise< boost::mpl::int_<0>() >
              > {};

          public:
            //! The result type.
            typedef boost::result_of< count(T_EXPR) > :: type type;

            //! The call itself.
            static type call( Sections<T_SQ> const& ) { return type; }
        };

      template<class T_EXPR> template<class T_SQ>
        struct Sections<T_EXPR> :: begin
        {
          //! The first iterator type.
          typedef section_iterator_begin< T_SQ > type;

          //! The call itself.
          static type call(T_SQ& sq) { return type(sq); }
        };

      template<class T_EXPR> template<class T_SQ>
        struct Sections<T_EXPR> :: end
        {
          //! The first iterator type.
          typedef section_iterator_end< T_SQ > type;

          //! The call itself.
          static type call(T_SQ& sq) { return type(sq); }
        };
    } // namespace sequences.
  } // namespace load_n_save
} // namespace LaDa


#endif
