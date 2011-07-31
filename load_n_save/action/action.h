#ifndef LADA_LNS_ACTION_H
#define LADA_LNS_ACTION_H

#include "LaDaConfig.h"

#include <boost/utility/enable_if.hpp>

#include "../string_type.h"

#include "action_base.h"
#include "type_to_regex.h"
#include "type_to_string.h"
#include "string_to_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace action_
    {
      namespace details
      {
         template<class T> struct dummy;
      }
      //! \brief Class for standard actions with default.
      //! \details This class will assign the value of an option to a a variable.
      //!          A string regex must be available through the external template
      //!          function LaDa::load_n_save::t_String
      //!          LaDa::load_n_save::TypeToRegex::apply() for the variable type.
      //!          Furthermore, a string must be convertible to the type through
      //!          the external template function bool
      //!          Lada::load_n_save::string_to_type(). Finally, the default
      //!          should be of a type convertible to the type of the variable.
      template< class T_TYPE, class T_DEFAULT >
        class Action : public ActionBase
        {
          BOOST_STATIC_ASSERT
          ((
            boost::is_convertible<T_DEFAULT, T_TYPE>::type::value 
          ));
          public:
            //! Type the action holds.
            typedef T_TYPE t_Type;
            //! Type of the default;
            typedef T_DEFAULT t_Default;
     
            //! Constructor.
            Action   ( t_Type &_var, T_DEFAULT const& _def )
                   : variable_(_var), default_(_def) {}
            //! Copy constructor.
            Action   ( Action const &_c )
                   : variable_(_c.variable_), default_(_c.default_) {}
            //! Virtual destructor.
            virtual ~Action() {};
            //! Parses a value to the action.
            virtual bool operator()( t_String const& _str ) const
              { return string_to_type_(_str, typename details::dummy<T_TYPE>::type() ); }
            //! Parses a value to the action.
            virtual t_String operator()() const
              { return TypeToRegex<t_Type>::apply(); }
            //! Puts default value into action.
            virtual bool assign_default() const { variable_ = default_; return true; }
            //! Prints to string.
            virtual t_String str() const { return TypeToString<t_Type>::apply(variable_); }
     
          private:
            bool string_to_type_(t_String const& _str, boost::mpl::bool_<false> const& ) const
              { return StringToType<t_Type>::apply(_str, variable_); }
            bool string_to_type_(t_String const& _str, boost::mpl::bool_<true> const& ) const
              { return StringToType<t_Type>::apply(_str, variable_); }
            //! Variable the action holds.
            t_Type &variable_;
            //! Variable the action holds.
            t_Default const& default_;
        };

      //! \brief Class for a standard action without default.
      //! \details This class will assign the value of an option to a a variable.
      //!          A string regex must be available through the external template
      //!          function LaDa::load_n_save::t_String
      //!          LaDa::load_n_save::TypeToRegex::apply() for the variable type.
      //!          Furthermore, a string must be convertible to the type through
      //!          the external template function bool
      //!          Lada::load_n_save::string_to_type().
      template< class T_TYPE >
        class Action<T_TYPE, boost::mpl::bool_<false> > : public ActionBase
        {
          public:
            //! Type the action holds.
            typedef T_TYPE t_Type;
     
            //! Constructor.
            Action   ( t_Type &_var, boost::mpl::bool_<false> const& )
                   : variable_(_var) {}
            //! Copy constructor.
            Action( Action const &_c ): variable_(_c.variable_) {}
            //! Virtual destructor.
            virtual ~Action() {};
            //! Parses a value to the action.
            virtual bool operator()( t_String const& _str ) const
              { return string_to_type_(_str, typename details::dummy<T_TYPE>::type() ); }
            //! Parses a value to the action.
            virtual t_String operator()() const 
              { return TypeToRegex<t_Type>::apply(); }
            //! Prints to string.
            virtual t_String str() const { return TypeToString<t_Type>::apply(variable_); }
     
          private:
            bool string_to_type_(t_String const& _str, boost::mpl::bool_<false> const& ) const
              { return StringToType<t_Type>::apply(_str, variable_); }
            bool string_to_type_(t_String const& _str, boost::mpl::bool_<true> const& ) const
              { return StringToType<t_Type>::apply(_str, variable_); }
            //! Variable the action holds.
            t_Type &variable_;
        };

      //! \brief Virtualizes calls to special actions with a default.
      //! \details Special actions are helper classes to perform complex actions
      //!          upon the value of a text option, such as options with finite
      //!          number of values, or pushing the value into a vector.
      template< class T_ACTION, class T_DEFAULT >
        struct SpecialAction : public ActionBase
        {
          public:
            //! Type of the default;
            typedef T_DEFAULT t_Default;
     
            //! Constructor.
            SpecialAction   ( T_ACTION const &_ac, T_DEFAULT const& _def )
                          : action_(_ac), default_(_def) {}
            //! Copy constructor.
            SpecialAction   ( SpecialAction const &_c )
                          : action_(_c.action_), default_(_c.default_) {}
            //! Virtual destructor.
            virtual ~SpecialAction() {};
            //! Parses a value to the action.
            virtual bool operator()( t_String const& _str ) const
              {  return action_(_str); }
            //! Parses a value to the action.
            virtual t_String operator()() const
              { return action_(); }
            //! Puts default value into action.
            virtual bool assign_default() const { return action_.assign_default(default_); }
            //! Puts default value into action.
            virtual t_String str() const { return action_.str(); }
     
          private:
            //! Variable the action holds.
            T_ACTION const &action_;
            //! Variable the action holds.
            t_Default const& default_;
        };

      //! \brief Virtualizes calls to special actions without a default.
      //! \details Special actions are helper classes to perform complex actions
      //!          upon the value of a text option, such as options with finite
      //!          number of values, or pushing the value into a vector.
      template< class T_ACTION >
        struct SpecialAction< T_ACTION, boost::mpl::bool_<false> >: public ActionBase
        {
          public:
            //! Type of the default;
            typedef boost::mpl::bool_<false> t_Default;
     
            //! Constructor.
            SpecialAction   ( T_ACTION const &_ac, t_Default const& )
                          : action_(_ac) {}
            //! Copy constructor.
            SpecialAction   ( SpecialAction const &_c )
                          : action_(_c.action_) {}
            //! Virtual destructor.
            virtual ~SpecialAction() {};
            //! Parses a value to the action.
            virtual bool operator()( t_String const& _str ) const
              {  return action_(_str); }
            //! Parses a value to the action.
            virtual t_String operator()() const
              { return action_(); }
            //! Puts default value into action.
            virtual t_String str() const { return action_.str(); }
     
          private:
            //! Variable the action holds.
            T_ACTION const &action_;
        };

      //! \brief Class for id-type actions.
      //! \details This action makes sure that an id is an id.
      class IdAction : public ActionBase
      {
        public:
          //! Type the action holds.
          typedef t_String t_Type;
     
          //! Constructor.
          IdAction(t_String const &_var) : variable_(_var) {}
          //! Copy constructor.
          IdAction(IdAction const &_c) : variable_(_c.variable_) {}
          //! Virtual destructor.
          virtual ~IdAction() {};
          //! Parses a value to the action.
          virtual bool operator()( t_String const& _str ) const { return true; }
          //! Parses a value to the action.
          virtual t_String operator()() const { return variable_; }
          //! Puts default value into action.
          virtual bool assign_default() const { return false; }
          //! Prints to string.
          virtual t_String str() const { return variable_; }
     
        private:
          //! Variable the action holds.
          t_String variable_;
      };

      namespace details
      {
         template<class T> struct dummy { typedef boost::mpl::bool_<false> type; };
      }
    } // namespace action_ 
  } // namespace load_n_save
} // namespace LaDa

#endif
