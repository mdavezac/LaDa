#ifndef LADA_MATH_SERIALIZE_H
#define LADA_MATH_SERIALIZE_H

#include "LaDaConfig.h"

#include <sstream>

#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/seq/elem.hpp>

#include "eigen.h"

#ifdef LADA_WITH_LNS
# include <load_n_save/action/string_to_type.h>
# include <load_n_save/action/type_to_string.h>
# include <load_n_save/action/type_to_regex.h>
#endif 

namespace boost {
  namespace serialization {

#   ifdef LADA_MACRO
#     error LADA_MACRO already defined.
#   endif
#   ifdef LADA_SEQX
#     error LADA_SEQX already defined.
#   endif
#   ifdef LADA_SEQY
#     error LADA_SEQY already defined.
#   endif
#   define LADA_SEQY (1)(2)(3)(4)(5)\
                     (1)(2)(3)(4)(5)\
                     (1)(2)(3)(4)(5)\
                     (1)(2)(3)(4)(5)
#   define LADA_SEQX (2)(2)(2)(2)(2)\
                     (3)(3)(3)(3)(3)\
                     (4)(4)(4)(4)(4)\
                     (5)(5)(5)(5)(5)
#   define LADA_MACRO(z, n, type) \
      template<class Archive>       \
        void serialize(Archive & _ar, \
            Eigen::Matrix<type, BOOST_PP_SEQ_ELEM(n, LADA_SEQX), BOOST_PP_SEQ_ELEM(n, LADA_SEQY)>& _g, \
            const unsigned int _version) \
        { \
          for(size_t i(0); i < BOOST_PP_SEQ_ELEM(n, LADA_SEQX); ++i) \
            for(size_t j(0); j < BOOST_PP_SEQ_ELEM(n, LADA_SEQY); ++j) \
              _ar & _g(i,j); \
        }
    BOOST_PP_REPEAT(20, LADA_MACRO, LaDa::types::t_real)
    BOOST_PP_REPEAT(20, LADA_MACRO, LaDa::types::t_int)
#   undef LADA_MACRO
#   undef LADA_SEQX
#   undef LADA_SEQY

    //! Serializes eigen real vectors.
    template<class Archive>
      void serialize(Archive & _ar, LaDa::math::Affine3d & _g, const unsigned int _version)
       { _ar & _g.matrix(); }
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
