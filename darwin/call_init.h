//
//  Version: $Id$
//
#ifndef _DARWIN_ALLOY_LAYERS_CALLINIT_H_
#define _DARWIN_ALLOY_LAYERS_CALLINIT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


namespace LaDa
{

  namespace GA
  {
    //! A functor which determines whether or not to call an init member function.
    //! \tparam T_POLICY policy for which to check wether is contains an
    //!                  init member function.
    //! \tparam T_CALLER type of the class containing the policy. This type
    //!                  is part of the init member function's signature.
    template< class T_POLICY, class T_CALLER >
      class CallInit
      {
        protected:
          //! \brief A metafunctor using SNIFAE to dectect whether a policy has
          //!        an init method.
          class HasInitMethod
          {
            //! \brief Structure can be instantiated only for certain class.
            //! \details Class U must have a member with the correct signature.
            //!          Function Test then determines the name of this member.
            template<class U, void (U::*)( const T_CALLER& )> struct SFINAE {};
            //! This overload can be instantiated only with the correct class and class member.
            template<class U> static char Test(SFINAE<U, &U::init >*);
            //! This overload can be instantiated with any class.
            template<class U> static int Test(...);
            public:
              //! The resulting value.
              static const bool value = sizeof(Test<T_POLICY>(0)) == sizeof(char);
          };
          //! A simple tag for overloading;
          template< bool docall > class Tag {};

        public:
          //! Type of the policy class.
          typedef T_POLICY t_Policy;
          //! Type of the caller class.
          typedef T_CALLER t_Caller;
          //! \brief Branches through an overloaded function.
          //! \details Casting to policy type is also handled at this point.
          static void call( t_Caller &_caller )
            { branch( _caller, _caller, Tag< HasInitMethod::value >() ); }

        protected:
          //! Does have init member.
          static void branch( t_Policy& _p,
                              const t_Caller & _caller,
                              const Tag< true > & ) { _p.init( _caller ); }
          //! Does not have init member.
          static void branch( t_Policy& _p,
                              const t_Caller & _caller,
                              const Tag< false > & ) { }
      };

  }
} // namespace LaDa

#endif
