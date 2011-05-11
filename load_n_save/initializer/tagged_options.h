//
//  Version: $Id: tagged_options.h 1227 2009-07-14 02:17:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_TAGGEDOPTIONS_H_
#define _LADA_LOADNSAVE_TAGGEDOPTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <boost/lambda/core.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/static_assert.hpp>

#include "../tree/tree.h"

#include "../transforms/gettag.h"
#include "../transforms/getname.h"
#include "../transforms/getoptions.h"
#include "../grammar/tags.h"
#include "../grammar/options.h"

#include "option_.h"

namespace LaDa 
{
  namespace load_n_save
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

    class TaggedOptions;

    //! Dumps id options to a stream.
    std::ostream& operator<<( std::ostream& _stream, TaggedOptions const &_idops );

  
    //! Gathers id-ed options and performs operations.
    class TaggedOptions
    {
      friend std::ostream& operator<<( std::ostream& _stream, TaggedOptions const &_idops ); 
      protected:
        //! Type of the container
        typedef std::vector< std::string > t_Options;

      public:
        //! Constructor.
        template< class T_EXPR >
          TaggedOptions   ( T_EXPR const &_expr,
                            size_t _tag = tags::option::default_ )
                        : options_( new t_Options ), tag_( _tag )
          {
            BOOST_STATIC_ASSERT
              (( boost::proto::matches<T_EXPR, grammar::Content>::type::value ));
            boost::fusion::for_each( transform::GetOptions()(_expr), fillvector(*options_, tag_) );
          }
        //! Constructor.
        TaggedOptions   ( boost::mpl::bool_<false> const&,
                          size_t _tag = tags::option::default_ )
                      : tag_( _tag ) {}
        //! Copy Constructor.
        TaggedOptions  ( TaggedOptions const &_c )
                       : tag_( _c.tag_ ), options_( _c.options_ ) {}

        //! Returns true if the expected options are in the tree.
        bool operator()( tree::Section const & _sec ) const;

        //! Returns the unknown sections in a tree.
        TaggedOptions unknown( tree::Section const & _sec ) const;

        //! return true if container is empty.
        bool empty() const { return options_ ? options_->size() == 0: true; }

      protected:
        //! Holds id options.
        boost::shared_ptr<t_Options> options_;
        //! Which tag we are trying to identify.
        size_t tag_;

        //! Fills a vector when option is identified.
        struct fillvector
        {
          fillvector( t_Options& _op, size_t _tag ) : vec(_op), tag(_tag) {}
          fillvector( fillvector& _c ) : vec(_c.vec), tag(_c.tag) {}
          typedef void result_type;
          template< class T_EXPR >
            void operator()( transform::OptionalGroups<T_EXPR> const &_expr ) const
              { boost::fusion::for_each( _expr.view, *this ); }
          template< class T_EXPR >
            typename boost::disable_if
            <
              typename details::is_view<T_EXPR> :: type,
              void
            > :: type operator()( T_EXPR const &_expr ) const
            {
              BOOST_STATIC_ASSERT
                (( boost::proto::matches<T_EXPR, grammar::OrnateOption>::type::value ));
              if( not tag )
                vec.push_back( transform::GetName()( _expr ) );
              else if( transform::GetTag()( _expr ) & tag )
                vec.push_back( transform::GetName()( _expr ) );
            }
          template< class T_EXPR >
            typename boost::enable_if
            <
              typename details::is_view<T_EXPR> :: type,
              void
            > :: type operator()( T_EXPR const &_expr ) const
            {
              if( vec.size() > 0 and vec.back() == ")" ) vec.back() = ") or ("; 
              else vec.push_back("(");
              boost::fusion::for_each( _expr, *this );
              if( vec.back() == "(" ) vec.pop_back();
              else if( vec.back() == ") or (" ) vec.back() = ")";
              else vec.push_back(")");
            }
          void operator()( boost::mpl::bool_<false> const & ) const {};
          t_Options &vec;
          size_t tag;
        };
    };

    template< class T_EXPR >
      bool id_options( tree::Section const & _sec, T_EXPR const& _content  )
      {
        BOOST_STATIC_ASSERT
        (( 
             boost::proto::matches<T_EXPR, grammar::OrnateOptions>::type::value
          or boost::proto::matches<T_EXPR, grammar::Content>::type::value
        ));
        namespace bl = boost::lambda;
        bool result(true);
        boost::fusion::for_each
        (
          transform::GetOptions()(_content),
          bl::var(result) &= bl::bind(initializer::Option_(_sec, true, tags::option::id), bl::_1)
        );
        return result;
      }

     inline bool id_options( tree::Section const &, boost::mpl::bool_<false> const& )
       { return true; }

  } // namespace load_n_save

} // namespace LaDa


#endif
