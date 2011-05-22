#ifndef LADA_LNS_ACTION_BASE_H
#define LADA_LNS_ACTION_BASE_H

#include "LaDaConfig.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace action_
    {
      //! Class to determine whether an acion is special or standard.
      template< class T_DATA >
        class is_special_action
        {
            template<class U, class T > struct SFINAE {};
            //! This overload can be instantiated only with the correct class and class member.
            template<class U> static char Test(SFINAE<U, typename U::action >*);
            //! This overload can be instantiated with any class.
            template<class U> static int Test(...);
          public:
            //! The resulting value.
            static const bool value = sizeof(Test<T_DATA>(0)) == sizeof(char);
        };
  
      //! Abstract base class of actions.
      struct ActionBase
      {
        //! Virtual destructor.
        virtual ~ActionBase() {};
        //! Parses a value to the action.
        virtual bool operator()( t_String const& ) const = 0;
        //! Parses a value to the action.
        virtual t_String operator()() const = 0;
        //! Assigns default value to action.
        virtual bool assign_default() const { return false; }
        //! Prints a value
        virtual t_String str() const = 0;
      };

    } // namespace action_.
                              
  } // namespace load_n_save
} // namespace LaDa

#endif
