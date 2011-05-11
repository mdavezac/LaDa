//
//  Version: $Id: options.h 1137 2009-05-22 23:35:00Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_OPTIONS_H_
#define _LADA_LOADNSAVE_TRANSFORMS_OPTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>
#include <boost/proto/tags.hpp>

#include <boost/fusion/include/single_view.hpp>
#include <boost/fusion/include/joint_view.hpp>
#include <boost/type_traits/add_const.hpp>

#include "../grammar/sections.h"
#include "terminals.h"
#include "false.h"
#include "protect.h"
#include "join.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;
      namespace fusion = boost::fusion;

      struct OptionsImpl;
      struct OptionWithAction;

      //! \details Parses a section.
      //! \brief this sidestep allows us to make sure a lonely
      //!        grammar::OrnateOptions is stuck within a vector of 1
      //!        component, rather then returned directly.
      struct Options : public proto::or_
        <
          proto::when< grammar::OrnateOption, protect( OptionWithAction() ) >,
          proto::otherwise< OptionsImpl() >
        >  {};

      //! Cases for the OptionsImpl.
      struct OptionsImplCases
      {
        //! Primary case matches nothing.
        template< class Tag > struct case_ : public proto::not_<proto::_> {};
      };

      //! Parses sections.
      struct OptionsImpl : public proto::switch_<OptionsImplCases> {};

      //! subscript (action) case.
      struct OptionWithAction : public proto::or_
        <
          proto::when
          < 
            proto::subscript
            < 
              grammar::option_terminal,
              proto::or_
              <
                proto::terminal<proto::_ const&>,
                proto::terminal<proto::_ &> 
              >
            >,
            join( OptionsImpl(proto::_left), false_(), OptionsImpl(proto::_right) )
          >,
          proto::when
          < 
            proto::subscript
            < 
              grammar::OrnateOption,
              proto::or_
              <
                proto::terminal<proto::_ const&>,
                proto::terminal<proto::_ &> 
              >
            >,
            join( OptionsImpl(proto::_left), OptionsImpl(proto::_right) )
          >,
          proto::when
          <
            proto::assign
            <
              grammar::option_terminal,
              grammar::Values
            >,
            join( OptionsImpl( proto::_left ), protect(OptionsImpl(proto::_right)), false_() )
          >,
          proto::otherwise
          <
            join( OptionsImpl(), false_(), false_() )
          >
        > {};

      //! terminal case.
      template<> struct OptionsImplCases::case_< proto::tag::terminal >
        : public proto::or_
        <
          proto::when
          < 
            grammar::option_terminal,
            protect( name(grammar::option_terminal) )
          >,
          proto::when
          < 
            proto::terminal<proto::_ const&>,
            protect_cref( call(proto::terminal<proto::_ const&>) )
          >,
          proto::when
          < 
            proto::terminal<proto::_ &>,
            protect_ref( var(proto::terminal<proto::_ &>) )
          >,
          proto::when
          < 
            proto::terminal<proto::_>,
            protect( call( proto::terminal<proto::_> ) )
          >
        > {};

      //! plus case.
      template<> struct OptionsImplCases::case_< proto::tag::plus >
        : public proto::or_
        <
          proto::when
          <
            proto::plus< grammar::Content, grammar::OrnateSection >,
            Options( proto::_left )
          >,
          proto::when
          <
            proto::plus< grammar::OrnateOption, grammar::OrnateOption >,
            join( Options(proto::_left), Options(proto::_right)  )
          >,
          proto::when
          <
            proto::plus< grammar::OrnateOptions, grammar::OrnateOption >,
            join( OptionsImpl(proto::_left), Options(proto::_right)  )
          >,
          proto::when
          <
            grammar::OrnateOptions,
            join( OptionsImpl(proto::_left), OptionsImpl(proto::_right)  )
          >,
          proto::when
          <
            proto::plus< grammar::Content, grammar::OrnateSection >,
            OptionsImpl( proto::_left )
          > 
        > {};

      //! assign case.
      template<> struct OptionsImplCases::case_< proto::tag::assign >
        : public proto::or_
        <
          proto::when
          <
            proto::assign
            < 
              grammar::option_terminal, 
              proto::terminal<proto::_> 
            >,
            join( OptionsImpl(proto::_left), protect(OptionsImpl(proto::_right))  )
          >,
          proto::when
          <
            proto::assign< grammar::option_terminal, grammar::Values >,
            join( OptionsImpl(proto::_left), protect(OptionsImpl(proto::_right))  )
          >
        > {};

      //! logical_or case.
      template<> struct OptionsImplCases::case_< proto::tag::logical_or >
        : public proto::when 
        <
          proto::logical_or< grammar::Values, grammar::Values >,
          join( OptionsImpl(proto::_left), OptionsImpl(proto::_right) )
        > {};

      //! == case.
      template<> struct OptionsImplCases::case_< proto::tag::equal_to >
        : public proto::when 
        <
          proto::equal_to< grammar::string_terminal, proto::_ >,
          protect( join( OptionsImpl(proto::_left), OptionsImpl(proto::_right) ) )
        > {};

    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif
