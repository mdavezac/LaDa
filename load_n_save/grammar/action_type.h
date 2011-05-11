//
//  Version: $Id: action_type.h 1231 2009-07-17 05:12:39Z davezac $
//

#ifndef _LADA_LOADNSAVE_GRAMMAR_ACTION_TYPE_H_
#define _LADA_LOADNSAVE_GRAMMAR_ACTION_TYPE_H_

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
        template< class T_TYPE > struct action_type {};
      }

      //! External type terminals.
      typedef proto::terminal< details::action_type<proto::_> > action_type;

    } // namespace grammar.


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
         boost::proto::term< details::action_type<T_ACTION> > \
       > 
  // include Dsel specialization.
  template<class T_ACTION>  
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
        typedef T_ACTION t_Type;
        //! The base type.
        typedef proto::extends
                < 
                  DSEL_TYPENAME,
                  Dsel< DSEL_TYPENAME >,
                  Domain 
                > t_Base;
        //! Holds a reference to the external type.
        T_ACTION var_;

      public:
        //! A general constructor to abstract out expression types.
        Dsel( T_ACTION const & _expr ) : t_Base(), var_(_expr){}
        //! Copy constructor.
        Dsel( Dsel const& _a ) : t_Base(_a), var_(_a.var_) {};
        //! \brief Serializes through \a _a.
        //! \details Dispatches call to either T_ACTION :: lns_access()
        //!           or load_n_save::lns_access(T_ARCHIVE const&, T_ACTION& )
        template< class T_ARCHIVE > 
          bool operator()( T_ARCHIVE const& _a ) const
            { return access<T_ARCHIVE, T_ACTION>()( _a, var_ ); }
    };

  namespace details
  {
    template< class T_ACTION >
      struct result_of_action
      {
        typedef Dsel< DSEL_TYPENAME > type;
      };
  }
  template< class T_ACTION >
    Dsel< DSEL_TYPENAME > const make_action( T_ACTION const &_type ) 
      { return Dsel< DSEL_TYPENAME >(_type); }
# undef DSEL_TYPENAME

#endif
