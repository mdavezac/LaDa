//
//  Version: $Id: initializer.h 1099 2009-05-10 06:31:25Z davezac $
//


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/xpressive/regex_primitives.hpp>

namespace LaDa 
{
  namespace load_n_save
  {
    namespace bx = boost::xpressive;
    
    const boost::xpressive::sregex XmlFormat::grammar_::option
      =    (bx::s1 =+bx::alnum) << -(*bx::_s) << '=' << -(*bx::_s)
        << (bx::s2=-!bx::as_xpr('\"'))
        << (bx::s3=+bx::alnum) << bx::s2;
    const boost::xpressive::sregex XmlFormat::grammar_::line_tag
      =    bx::as_xpr('<') << (bx::s1=+bx::alnum)
        << *( options || bx::_ln || bx::_s )
        << bx::as_xpr("/>") 
    const boost::xpressive::sregex XmlFormat::grammar_::start_tag
      =    bx::as_xpr('<') << (bx::s1=+bx::alnum)
        << *( options || bx::_ln || bx::_s )
        << bx::as_xpr('>') 
    const boost::xpressive::sregex XmlFormat::grammar_::end_tag
      =    bx::as_xpr("/") << (bx::s1=+bx::alnum)
        << *bx::_s << bx::as_xpr('>') 
    const boost::xpressive::sregex XmlFormat::grammar_::tag;
      =    ( line_tag || ( start_tag << *tag << end_tag ) );

    boost::xpressive::sregex grammar_::option_regex( const std::string& _name )
    {
      return    bx::as_xpr(_name) << -(*bx::_s) << '=' << -(*bx::_s)
             << (bx::s2=-!bx::as_xpr('\"'))
             << (bx::s3=+bx::alnum) << bx::s2;
    }

    boost::xpressive::sregex grammar_::option_regex( const std::string& _name,
                                                     const std::string& _value)
    {
      return    bx::as_xpr(_name) << -(*bx::_s) << '=' << -(*bx::_s)
             << (bx::s2=-!bx::as_xpr('\"'))
             << bx::as_xpr(value) << bx::s2;
    }
  } // namespace load_n_save
} // namespace LaDa

