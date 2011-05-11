//
//  Version: $Id: grammar.h 1215 2009-07-08 00:28:08Z davezac $
//

#ifndef _LADA_LOADNSAVE_GRAMMAR_H_
#define _LADA_LOADNSAVE_GRAMMAR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _LADADEBUG
# include <iostream>
#endif
#include <boost/fusion/include/at_c.hpp>

#include <boost/proto/tags.hpp>
#include <boost/proto/operators.hpp>
#include <boost/proto/domain.hpp>

#include "dsel.h"
#include "terminals.h"
#include "options.h"
#include "sections.h"

#include <opt/types.h>

namespace LaDa 
{
  namespace load_n_save
  {
    namespace grammar
    {
      namespace proto = boost::proto;
     
      //! Cases for the grammar.
      struct GrammarCases
      {
          //! The primary template matches nothing.
          template<typename Tag> struct case_ : proto::not_<proto::_> {};
      };
     
      //! The allowable grammars.
      struct Grammar : public proto::switch_<GrammarCases> {};
     
      //! Grammar case: terminals.
      template<> struct GrammarCases::case_<proto::tag::terminal> :
        public proto::terminal< proto::_ > {};
      //! Grammar case: plus.
      template<> struct GrammarCases::case_<proto::tag::plus> :
        public proto::or_
        <
          proto::plus< OrnateOptions, OrnateOptions >,
          proto::plus< OrnateOptions, OrnateSections >,
          proto::plus< Content, OrnateSections >,
          proto::plus< OrnateSections, OrnateSections >
        > {};
    
      //! Grammar case: terminals.
      template<> struct GrammarCases::case_<proto::tag::assign> : 
        public proto::or_
        <
          proto::assign< option_function, TypedValues<proto::_> >,
          proto::assign< section_function, Content >,
          tag_grammar,
          help_grammar,
          default_grammar,
          action_grammar
        > {};
     
      //! Grammar case: ==.
      template<> struct GrammarCases::case_<proto::tag::equal_to> : public Value {};
     
      //! Grammar case: ||.
      template<> struct GrammarCases::case_<proto::tag::logical_or> : public proto::or_
        <
          Values,
          proto::logical_or< OrnateOptions, OrnateOptions >,
          proto::logical_or< OrnateSections, OrnateSections >
        > {};
     
      //! Grammar case: subscript.
      template<> struct GrammarCases::case_<proto::tag::function> :
        public proto::or_
        <
          option_function,
          section_function
        > {};
     
      //! Domain of the grammar.
      struct Domain : public proto::domain< proto::pod_generator< Dsel >, Grammar > {};
    } // namespace grammar
    

  } // namespace load_n_save

} // namespace LaDa


#endif
