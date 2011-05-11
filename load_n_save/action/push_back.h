//
//  Version: $Id: push_back.h 1250 2009-07-26 21:04:07Z davezac $
//

#ifndef _LADA_LOADNSAVE_PUSH_BACK_ACTION_H_
#define _LADA_LOADNSAVE_PUSH_BACK_ACTION_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include "grammar/action_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    //! Pushes back into a vector.
    template<class T_CONTAINER>
      class PushBack
      {
        public:
          //! A trait wich says this is an action.
          typedef void action;
          //! The return type is bool.
          typedef bool result_type;
          //! Constructor. 
          PushBack(T_CONTAINER &_cont) : container_(_cont) {}
          //! Copy
          PushBack(PushBack const &_c) : container_(_c.container_) {}
          //! Does bitwise operation.
          result_type operator()( std::string const &_x ) const
          { 
            typename T_CONTAINER::value_type v;
            if( not initializer::parse_value( _x, v ) ) return false;
            container_.push_back( v );
            return true; 
          }

          template< class T_ARCHIVE >
            bool lns_access( T_ARCHIVE &_ar )
            {
              size_t i(0);
              typename T_CONTAINER::value_type v;
              for(; not (_ar & grammar::ext( v ) ); ++i )  container_.push_back( v );
              return i > 0;
            } 

        protected:
          //! Holds reference to container.
          T_CONTAINER &container_;
      };

    //! Pushes back into a vector.
    template<class T_CONTAINER, class T_CONVERTER>
      class PushBackWithConverter : PushBack<T_CONTAINER>
      {
        public:
          //! A trait wich says this is an action.
          typedef void action;
          //! The return type is bool.
          typedef bool result_type;
          //! Constructor. 
          PushBackWithConverter   ( T_CONTAINER &_cont, T_CONVERTER &_converter )
                                : PushBack<T_CONTAINER>(_cont), converter_(_converter) {} 
          //! CopyConstructor. 
          PushBackWithConverter   ( PushBackWithConverter const &_c )
                                : PushBack<T_CONTAINER>(_c.container_), converter_(_c.converter_) {} 
          //! Does bitwise operation.
          result_type operator()( std::string const &_x ) const
          { 
            if( not initializer::parse_value( _x, converter_ ) ) return false;
            PushBack<T_CONTAINER>::container_.push_back( converter_ );
            return true; 
          }

          template< class T_ARCHIVE >
            bool lns_access( T_ARCHIVE &_ar )
            {
              size_t i(0);
              for(; not (_ar & grammar::ext( converter_ ) ); ++i )
                PushBack<T_CONTAINER>::container_.push_back( converter_ );
              return i>0;
            } 

        protected:
          //! Holds reference to an action which converts to the container value type.
          T_CONVERTER &converter_;
      };

    //! Returns an action to push_back items into a container.
    template< class T_CONTAINER >
      PushBack<T_CONTAINER> push_back(T_CONTAINER &_cont)
        { return PushBack<T_CONTAINER>(_cont); }

#   ifdef DSEL_TYPENAME
#     error DSEL_TYPENAME already exists.
#   endif
#   ifdef T_TYPE
#     error T_TYPE already exists
#   endif
#   define T_TYPE 
#   define DSEL_TYPENAME                                         \
         boost::proto::expr                                      \
         <                                                       \
           boost::proto::tag::terminal,                          \
           boost::proto::term< grammar::details::external_type<const T_TYPE> >  \
         > 
    //! Returns an action to push_back items into a container.
    template< class T_CONTAINER, class T_CONVERTER >
      typename grammar::details::result_of_action< PushBackWithConverter<T_CONTAINER, T_CONVERTER> > :: type
        push_back( T_CONTAINER &_cont, T_CONVERTER& _conv )
          { return grammar::make_action( PushBackWithConverter<T_CONTAINER, T_CONVERTER>(_cont, _conv) ); }

#   undef T_TYPE
#   undef DSEL_TYPENAME   
  }
}

#endif
