//
//  Version: $Id: transforms.h 1101 2009-05-10 22:29:36Z davezac $
//

#ifndef _LADA_LOADNSAVE_TRANSFORMS_H_
#define _LADA_LOADNSAVE_TRANSFORMS_H_

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

#include "terminals.h"
#include "push_back_section.h"
#include "section_sum.h"

#ifdef _LADADEBUG
  struct display
  {
      template<typename T>
      void operator()(T const &t) const; 
  };
  template< typename T > void display :: operator()( T const &t ) const 
  {  boost::proto::display_expr( t ); }
  template<> void display :: operator()( boost::fusion::nil const & ) const 
  {  std::cout << " nil\n"; }
#endif

namespace LaDa 
{
  namespace load_n_save
  {
    namespace transform
    {
      namespace proto = boost::proto;
      namespace fusion = boost::fusion;

      //! Cases for the transform.
      struct transform_cases
      {
        //! Primary case matches nothing.
        template< class Tag > struct case_ : public proto::not_<proto::_> {};
      };

      struct transform : public proto::switch_<transform_cases> {};

      //! terminal case.
      template<> struct transform_cases::case_< proto::tag::terminal >
        : public proto::or_
        <
          proto::when< section_terminal, name( section_terminal ) >,
          proto::when< option_terminal, name(section_terminal) >
        > {};

      //! assign case.
      template<> struct transform_cases::case_< proto::tag::assign >
        : public proto::when
        <
          proto::assign< section_terminal, option_terminal >,
          fusion::joint_view
          <
            boost::add_const< transform( proto::_left ) >,
            boost::add_const< transform( proto::_right ) >
          >( transform( proto::_left ), transform( proto::_right ) )
        > {};

      //! assign case.
      template<> struct transform_cases::case_< proto::tag::plus >
        : public proto::or_
        <
          proto::when
          <
            proto::plus< OrnateSection, OrnateSection >,
            section_sum( transform(proto::_left), transform(proto::_right) )
          >,
          proto::when
          <
            proto::plus< OrnateSections, OrnateSection >,
            push_back_section( transform(proto::_left), transform(proto::_right) )
          >,
          proto::when
          <
            OrnateSections,
            fusion::joint_view
            < 
              boost::add_const<transform( proto::_left )>,
              boost::add_const<transform( proto::_right )>
            >( transform(proto::_left), transform(proto::_right)  )
          >
        > {};
    }

 
  } // namespace load_n_save

} // namespace LaDa


#endif
