//
//  Version: $Id: external_type.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_GRAMMAR_EXTERNAL_TYPE_H_
#define _LADA_LOADNSAVE_GRAMMAR_EXTERNAL_TYPE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _LADADEBUG
# include <iostream>
#endif
#include <boost/fusion/include/at_c.hpp>

#include <boost/proto/tags.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/equal_to.hpp>

#include <opt/types.h>

#include "dsel.h"

namespace LaDa 
{
  namespace load_n_save
  {

    namespace grammar
    {
      namespace proto = boost::proto;

      namespace details
      {
        //! Tag for serialization for external types.
        template< class T_TYPE > struct external_type {};

        template< class T_TYPE >
          struct dsel_ext_type : public proto::terminal
                <
                  external_type<T_TYPE>
                > :: type {};

      }

      //! External type terminals.
      typedef proto::terminal< details::external_type<proto::_> > external_type;


    } // namespace grammar.

    //! \brief Access granted to archives.
    //! \details boost::serialization has a better way of doing this.
    //!          The problem is that if we declare a template function here, 
    //!          The look up will happen here as well, rather than wait for
    //!          instantiation. Which means the template function cannot be
    //!          effectively overloaded. The solution relies on some subtle
    //!          overloading of the function for the third argument giving the
    //!          archive version. However, we have no third argument here. This
    //!          leads to requiring overly generous access as a friend
    //!          template class. And/or to a bit more boiler plate.
    //! \todo Find some way to use a function rather than a class, like
    //!       boost::serialization.
    template<class T_ARCHIVE, class T_TYPE>
      struct access
      {
        bool operator()( T_ARCHIVE const& _ar, T_TYPE &_atom ) const
        {
          _atom.lns_access( _ar ); 
        }
      };


  } // namespace load_n_save
} // namespace LaDa


#elif defined( FROM_LADA_DSEL )

// g++ bug?
// It seems that using a typedef in templated class as a partial specialization
// will not compile. However, using the type directly will.
# ifdef DSEL_TYPENAME
#   error DSEL_TYPENAME already exists.
# endif
# define DSEL_TYPENAME                                         \
       boost::proto::expr                                      \
       <                                                       \
         boost::proto::tag::terminal,                          \
         boost::proto::term< details::external_type<T_TYPE> > \
       > 
  // include Dsel specialization.
  template<class T_TYPE>  
    struct Dsel< DSEL_TYPENAME >
      : public proto::extends
               < 
                 DSEL_TYPENAME,
                 Dsel< DSEL_TYPENAME >,
                 Domain 
               >
    {
      protected:
        //! External type.
        typedef T_TYPE t_Type;
        //! The base type.
        typedef proto::extends
                < 
                  DSEL_TYPENAME,
                  Dsel< DSEL_TYPENAME >,
                  Domain 
                > t_Base;
        //! Holds a reference to the external type.
        T_TYPE &var_;

      public:
        //! A general constructor to abstract out expression types.
        Dsel( T_TYPE & _expr ) : t_Base(), var_(_expr){}
        //! Copy constructor.
        Dsel( Dsel const& _a ) : t_Base(_a), var_(_a.var_) {};
        //! \brief Serializes through \a _a.
        //! \details Dispatches call to either T_TYPE :: lns_access()
        //!           or load_n_save::lns_access(T_ARCHIVE const&, T_TYPE& )
        template< class T_ARCHIVE > 
          bool operator()( T_ARCHIVE const& _a ) const
            { return access<T_ARCHIVE, T_TYPE>()( _a, var_ ); }
    };
  template< class T_TYPE >
    Dsel< DSEL_TYPENAME > const ext( T_TYPE &_type ) 
      { return Dsel< DSEL_TYPENAME >(_type); }
# undef DSEL_TYPENAME

#endif
