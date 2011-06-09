#ifndef LADA_MATH_SERIALIZE_H
#define LADA_MATH_SERIALIZE_H

#include "LaDaConfig.h"

#include <sstream>

#include "eigen.h"

#ifdef LADA_WITH_LNS
# include <load_n_save/action/string_to_type.h>
# include <load_n_save/action/type_to_string.h>
# include <load_n_save/action/type_to_regex.h>
#endif 

namespace boost {
  namespace serialization {

    //! Serializes eigen real vectors.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::rVector3d & g, const unsigned int version)
     { ar & g.x(); ar & g.y(); ar & g.z(); }
     //! Serializes eigen integer vectors.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::iVector3d & g, const unsigned int version)
     { ar & g.x(); ar & g.y(); ar & g.z(); }
    //! Serializes eigen real matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::rMatrix3d & g, const unsigned int version)
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          ar & g(i,j);
    }
    //! Serializes eigen integer matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::iMatrix3d & g, const unsigned int version)
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          ar & g(i,j);
    }
  }
}

#ifdef LADA_WITH_LNS
  namespace LaDa
  {
    namespace load_n_save
    {
      //! Regex for vectors.
      template<> struct TypeToRegex<math::rVector3d, void>
      {
        //! Returns regex string.
        static t_String apply() 
          {
            return   TypeToRegex<types::t_real>::apply() + "\\s+"
                   + TypeToRegex<types::t_real>::apply() + "\\s+"
                   + TypeToRegex<types::t_real>::apply();
          }
      };

      //! Parses an math::rVector3d.
      template<> struct StringToType<math::rVector3d, void>
      {
        //! Functor.
        static bool apply( t_String const& _string, math::rVector3d &_value )
        {
          namespace bx = boost::xpressive;
          bx::smatch what;
          bx::sregex const sig = (bx::set= '+', '-');
          bx::sregex const exp =    bx::as_xpr('.')
                                 >> !( (bx::set='d','D','e','E') >> !sig)
                                 >> +bx::_d;
          bx::sregex rep = (bx::set='d','D','E');  \
          bx::sregex re =    (bx::s1 = (!sig >> +bx::_d >> !exp) ) >> +bx::_s 
                          >> (bx::s2 = (!sig >> +bx::_d >> !exp) ) >> +bx::_s 
                          >> (bx::s3 = (!sig >> +bx::_d >> !exp) );
          if( not bx::regex_match( _string, what, re ) ) return false;
          std::string const a( bx::regex_replace( std::string(what[1]), rep, std::string("e") ) );
          std::string const b( bx::regex_replace( std::string(what[2]), rep, std::string("e") ) );
          std::string const c( bx::regex_replace( std::string(what[3]), rep, std::string("e") ) );
          _value[0] = boost::lexical_cast<types::t_real>( a );
          _value[1] = boost::lexical_cast<types::t_real>( b );
          _value[2] = boost::lexical_cast<types::t_real>( c );

          return true;
        }
      };
      template<> struct TypeToString<math::rVector3d, void>
      {
        static t_String apply(math::rVector3d const &_value)
        {
          std::ostringstream sstr;
          sstr << _value(0) << " " << _value(1) << " " << _value(2);
          return sstr.str();
        }
      };
    }
  }
#endif

#endif
