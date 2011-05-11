//
//  Version: $Id: access.h 1250 2009-07-26 21:04:07Z davezac $
//

#ifndef _LADA_LNS_ACCESS_H_
#define _LADA_LNS_ACCESS_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace LaDa 
{
  namespace load_n_save
  {
    namespace parser_base
    {
      struct Section;
    }
    //! Non-intrusive external access. 
    template<class T_DATA> struct lns_access;

    //! \brief A metafunctor using SNIFAE to dectect whether a type has
    //!        an lns_access method.
    template< class T_DATA >
      class Access
      {
          //! \brief Structure can be instantiated only for certain class.
          //! \details Class U must have a member with the correct signature.
          //!          Function Test then determines the name of this member.
          template<class U, bool (U::*)( parser_base::Section const & )> struct SFINAE {};
          //! This overload can be instantiated only with the correct class and class member.
          template<class U> static char Test(SFINAE<U, &U::lns_access >*);
          //! This overload can be instantiated with any class.
          template<class U> static int Test(...);
          //! The resulting value.
          static const bool value = sizeof(Test<T_DATA>(0)) == sizeof(char);
        public:
          //! Makes member or external function call.
          static bool call( parser_base::Section const& _parser, T_DATA& _data ) 
            { return call_<T_DATA>( _parser, _data, boost::mpl::bool_<value>() ); }
        protected:
          //! Calls member function.
          template<class T> 
            static bool call_( parser_base::Section const& _parser, 
                               T& _data, boost::mpl::bool_<true> ) 
             { return _data.lns_access( _parser ); }
          //! External function.
          template< class T> 
            static bool call_( parser_base::Section const& _parser,
                               T& _data, boost::mpl::bool_<false> ) 
             { return lns_access<T_DATA>()( _parser, _data ); }
      };
  } // namespace load_n_save
} // namespace LaDa

#endif
