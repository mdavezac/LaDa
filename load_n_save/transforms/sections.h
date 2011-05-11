//
//  Version: $Id: sections.h 1170 2009-06-09 04:17:08Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_SECTIONS_H_
#define _LADA_LOADNSAVE_TRANSFORMS_SECTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/proto/core.hpp>
#include <boost/proto/transform.hpp>
#include <boost/proto/transform/when.hpp>
#include <boost/proto/matches.hpp>
#include <boost/proto/tags.hpp>

#include <boost/fusion/include/single_view.hpp>
#include <boost/fusion/include/joint_view.hpp>
#include <boost/type_traits/add_const.hpp>

#include "../grammar/sections.h"
#include "terminals.h"
#include "protect.h"
#include "join.h"
#include "options.h"
#include "false.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;
      namespace fusion = boost::fusion;

      struct SectionsVector : public proto::or_
        <
          proto::when< grammar::OrnateSection, protect( proto::_ ) >,
          proto::when
          <
            proto::plus< grammar::OrnateOptions, grammar::OrnateSections >,
            SectionsVector(proto::_right) 
          >,
          proto::when
          <
            proto::plus< grammar::Content, grammar::OrnateSections >,
            join( SectionsVector(proto::_left), SectionsVector(proto::_right)  )
          >,
          proto::when
          <
            proto::plus< grammar::OrnateSections, grammar::OrnateSections >,
            join( SectionsVector(proto::_left), SectionsVector(proto::_right) )
          >
        > {};


      struct SectionsImpl;
      struct SectionWithAction;
      struct Content;

      //! \details Parses a section.
      //! \brief this sidestep allows us to make sure lonely
      //!        grammar::OrnateSections are stuck within a vector of 1
      //!        component, rather then returned directly.
      struct Sections : public proto::or_
        <
          proto::when
          <
            grammar::OrnateSection,
            protect( SectionWithAction() )
          >,
          proto::otherwise< SectionsImpl() >
        >  {};
     
      //! Cases for the Sections.
      struct SectionsImplCases
      {
        //! Primary case matches nothing.
        template< class Tag > struct case_ : public proto::not_<proto::_> {};
      };

      //! Implementation of section parsing.
      struct SectionsImpl : public proto::switch_<SectionsImplCases> {};

      //! subscript (action) case.
      struct SectionWithAction : public proto::or_
        <
          proto::when 
          < // naked section with action.
            proto::subscript
            < 
              grammar::section_terminal,
              proto::terminal<proto::_ const&> 
            >,
            join
            ( 
              SectionsImpl(proto::_left),
              false_(),
              false_(),
              SectionsImpl(proto::_right) 
            )
          >,
          proto::when
          < // ornate section with action
            proto::subscript
            < 
              grammar::OrnateSection,
              proto::terminal<proto::_ const&> 
            >,
            join( SectionsImpl(proto::_left), SectionsImpl(proto::_right) )
          >,
          proto::when
          < // naked section no actions
            grammar::section_terminal,
            join( SectionsImpl(), false_(), false_(), false_() )
          >,
          proto::otherwise< join( SectionsImpl(), false_() ) >
        > {};

      struct Content : public proto::or_
        <
          proto::when // options alone.
          <
            grammar::OrnateOptions,
            boost::fusion::joint_view
            <
              boost::add_const
              <
                boost::fusion::single_view< boost::add_const<OptionsImpl()> >
              >,
              boost::add_const< false_() >
            >( boost::fusion::single_view< boost::add_const<OptionsImpl()> >(OptionsImpl()),
               false_() )
          >
        > {};

      //! terminal case.
      template<> struct SectionsImplCases::case_< proto::tag::terminal >
        : public proto::or_
        <
          proto::when
          < 
            grammar::section_terminal,
            protect( name(grammar::section_terminal) )
          >,
          proto::when
          < 
            proto::terminal<proto::_ const&>,
            protect_cref( call(proto::terminal<proto::_ const&>) )
          >
        > {};

      //! sum case.
      template<> struct SectionsImplCases::case_< proto::tag::plus >
        : public proto::or_
        <
          proto::when
          <
            proto::plus< grammar::OrnateSection, grammar::OrnateSection >,
            join( Sections(proto::_left), Sections(proto::_right) )
          >,
          proto::when
          <
            proto::plus< grammar::OrnateSections, grammar::OrnateSection >,
            join( SectionsImpl(proto::_left), Sections(proto::_right)  )
          >,
          proto::when
          <
            grammar::OrnateSections,
            join( SectionsImpl(proto::_left), SectionsImpl(proto::_right)  )
          > ,
          proto::when
          <
            proto::plus< grammar::OrnateOptions, grammar::OrnateSection >,
            Sections( proto::_right )
          >,
          proto::when 
          <
            proto::plus< grammar::Content, grammar::OrnateSection >,
            join( SectionsImpl(proto::_left), Sections(proto::_right) )
          > 
        > {};
  
      //! assign case.
      template<> struct SectionsImplCases::case_< proto::tag::assign >
        : public proto::or_
        <
          proto::when // section = option + option + ..  no subsections.
          <
            proto::assign< grammar::section_terminal, grammar::OrnateOptions >,
            join
            ( 
              SectionsImpl(proto::_left), 
              protect( Options(proto::_right) ), 
              false_() 
            )
          >,
          proto::when // section = section + section + ... no options.
          <
            proto::assign< grammar::section_terminal, grammar::OrnateSections >,
            join
            ( 
              SectionsImpl(proto::_left), 
              false_(),
              protect( Sections(proto::_right) )
            )
          >,
          proto::when
          < // in this case, grammar::Content means both options + sections.
            proto::assign< grammar::section_terminal, grammar::Content >,
            join
            (
              SectionsImpl(proto::_left),
              protect( Options(proto::_right) ),
              protect( Sections(proto::_right) )
            )
          >
        > {};
  
  
    } // namespace transform
 
  } // namespace load_n_save

} // namespace LaDa


#endif
