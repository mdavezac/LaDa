//
//  Version: $Id: tagged_sections.h 1227 2009-07-14 02:17:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_TAGGEDSECTIONS_H_
#define _LADA_LOADNSAVE_TAGGEDSECTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/static_assert.hpp>

#include "../tree/tree.h"

#include "../transforms/gettag.h"
#include "../transforms/getname.h"
#include "../transforms/getsections.h"
#include "../grammar/tags.h"
#include "../grammar/sections.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace fusion = boost::fusion;

    class TaggedSections;

    //! Dumps id options to a stream.
    std::ostream& operator<<( std::ostream& _stream, TaggedSections const &_idops );

  
    //! Gathers id-ed options and performs operations.
    class TaggedSections
    {
      friend std::ostream& operator<<( std::ostream& _stream, TaggedSections const &_idops ); 
      protected:
        //! Type of the container
        typedef std::vector< std::string > t_Sections;

      public:
        //! Constructor.
        template< class T_EXPR >
          TaggedSections   ( T_EXPR const &_expr, size_t _tag = tags::section::default_ ) 
                         : sections_( new t_Sections ), tag_( _tag )
          {
            BOOST_STATIC_ASSERT
              (( boost::proto::matches<T_EXPR, grammar::Content>::type::value ));
            boost::fusion::for_each( transform::GetSections()(_expr), fillvector(*sections_, tag_) );
          }
        //! Constructor.
        TaggedSections  ( boost::mpl::bool_<false> const&, size_t _tag = tags::section::default_ )
                       : tag_( _tag ) {}
        //! Copy Constructor.
        TaggedSections  ( TaggedSections const &_c )
                       : tag_( _c.tag_ ), sections_( _c.sections_ ) {}

        //! Checks that a tree has the sections stored in sections_
        bool operator()( tree::Section const & _sec ) const;

        //! Returns the unknown sections in a tree.
        TaggedSections unknown( tree::Section const & _sec ) const;
        
        //! return true if container is empty.
        bool empty() const { return sections_ ? sections_->size() == 0: true; }


      protected:
        //! Holds id options.
        boost::shared_ptr<t_Sections> sections_;
        //! Which tag we are trying to identify.
        size_t tag_;

        //! Fills a vector when section is identified.
        struct fillvector
        {
          fillvector( t_Sections& _op, size_t _tag ) : vec(_op), tag(_tag) {}
          fillvector( fillvector& _c ) : vec(_c.vec), tag(_c.tag) {}
          template< class T_EXPR >
            void operator()( T_EXPR const &_expr ) const
            {
              BOOST_STATIC_ASSERT
                (( boost::proto::matches<T_EXPR, grammar::OrnateSection>::type::value ));
              if( not tag )
                vec.push_back( transform::GetName()( _expr ) );
              else if( transform::GetTag()( _expr ) & tag )
                vec.push_back( transform::GetName()( _expr ) );
            }
          void operator()( boost::mpl::bool_<false> const & ) const {};
          t_Sections &vec;
          size_t tag;
        };
    };

  } // namespace load_n_save

} // namespace LaDa


#endif
