#include "LaDaConfig.h"

#include  <boost/xpressive/regex_algorithms.hpp>
#include  <boost/xpressive/regex_compiler.hpp>
#include  <boost/algorithm/string/trim.hpp>

#include "../tags.h"
#include "../exceptions.h"
#include "load.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace load
    {

      bool Load::Option:: operator()(t_String const& _name, xpr::Option const& _op ) const
      {
        tree::Section::const_iterator::option i_first = tree_.options_begin(_op.name);
        tree::Section::const_iterator::option const i_end = tree_.options_end(_op.name);

        bool found(false);
        bool parse_error(false);
        for(; i_first != i_end; ++i_first)
        {
          if( not i_first->fresh ) continue;
          found = true;

          // attempts parsing.
          parse_error = not parse_( _op.action(), i_first->value );
          if( grammar_only_ ) break;

          // does parsing.
          i_first->fresh = false;
          if( not parse_error ) _op.action( i_first->value );
          break;
        } 

        bool const is_id( _op.tag & load_n_save::idoption );
        bool const is_required( _op.tag & load_n_save::required );
        if( found )
        {
          if(parse_error and (not is_id))
            BOOST_THROW_EXCEPTION
            ( 
               error::option_parse_error() << error::option_name(_op.name)
                                           << error::section_name(_name) 
            );
          return not parse_error;
        }
        if( is_required )
          BOOST_THROW_EXCEPTION
          ( 
             error::required_option_not_found() << error::option_name(_op.name)
                                                << error::section_name(_name) 
          );
//       else if(is_required) return false;
        if( is_id ) return false;

        // assigns default value if it exits.
        if( not grammar_only_ ) _op.assign_default(); 

        return true;
      }

      bool Load::Option:: parse_( t_String const& _regex, t_String const &_value ) const
      {
        namespace bx = boost::xpressive;
        namespace ba = boost::algorithm;
        return bx::regex_match( ba::trim_copy(_value), bx::sregex::compile(ba::trim_copy(_regex)) );
      }

    } // namespace initializer.
  } // namespace load_n_save
} // namespace LaDa

