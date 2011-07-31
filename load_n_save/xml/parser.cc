#include "LaDaConfig.h"

#include <iostream>

#include <boost/lambda/lambda.hpp>
#include <boost/xpressive/regex_compiler.hpp>
#include <boost/xpressive/match_results.hpp>
#include <boost/xpressive/regex_primitives.hpp>
#include <boost/xpressive/regex_algorithms.hpp>

#include "../tree/tree.h"
#include "../exceptions.h"

#include "parser.h"
#include "tags.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xml
    {
 
      bool parse_options( tree::Section &_section,
                          std::string::const_iterator _first,  
                          std::string::const_iterator const _last,
                          size_t _line, tags::read const _tags, 
                          std::string const & _tag_name )
      {
        if( _first == _last ) return true;
        std::string::const_iterator const ofirst( _first );
        namespace bx = boost::xpressive; 
        bx::sregex const string_regex = bx::as_xpr('\"') >> *bx::_s >> (bx::s1=*(~bx::as_xpr('\"')))
                                                         >> *bx::_s >> '\"';
        bx::sregex const option 
          = ( (bx::s1=+bx::alnum) >> *bx::_s >> !( '=' >> *bx::_s >> (bx::s2=(string_regex | +(bx::alnum | '.'))) ) );
        bx::sregex const empty = *(bx::_s | bx::_n);
        
        bx::smatch match;
        bool result = true;
            
 
        do
        {
          if( bx::regex_match( _first, _last, empty ) ) break;
          bool const matched( bx::regex_search( _first, _last, match, option ) );
          if( not matched )
          { 
            if( _tags & tags::VALIDATE )
            {
              std::ostringstream sstr;
              size_t const l( _line + std::count( ofirst, _first, '\n' ) );
              std::for_each( _first, _last, sstr << boost::lambda::_1 );
              BOOST_THROW_EXCEPTION
              ( 
                 error::unparsable_option() << error::section_name(_tag_name) 
                                            << error::option_name(sstr.str())
                                            << error::line_number(l) 
              );
            }
            result = false;
            break;
          }
          if(     (_tags & tags::VALIDATE)
              and (not bx::regex_match( match.prefix().first, match.prefix().second, empty)) )
          {
            size_t const l( _line + std::count( ofirst, match.prefix().first, '\n' ) );
            BOOST_THROW_EXCEPTION
            ( 
               error::no_end_tag() << error::section_name(_tag_name) 
                                   << error::line_number(l) 
            );
            result = false;
          }
 
          _first = match[0].second;
          if( _tags & tags::READ )
          {
            std::string const name( match[1].str() );
            std::string const value( match[2].matched ? match[2].str(): "" );
            if( bx::regex_match( value, match, string_regex ) )
              _section.push_back( name, match[1].str() );
            else _section.push_back( name, value ); 
          }
        } while( _first != _last );
        return result;
      }

      bool parse( tree::Base &_base,
                  std::string::const_iterator _first,
                  std::string::const_iterator const _last,
                  tags::read const _tags, size_t _line )
      {
        std::string::const_iterator const ofirst( _first );
        namespace bx = boost::xpressive; 
        bx::sregex const start_tag = '<' >> (bx::s1=+bx::alnum);
        bx::sregex const empty = *(bx::_s | bx::_n);
        bx::smatch match;
        bool result = true;
     
        // on first entering loop, line must be empty,
        // then it must contain the last line looked at.
        while( _first != _last )
        {
          if( bx::regex_match( _first, _last, empty ) ) break;
          if( not bx::regex_search( _first, _last, match, start_tag ) )
          {
            if( not ( _tags & tags::VALIDATE ) ) break;
            size_t const l( _line + std::count( ofirst, _first, '\n' ) );
            BOOST_THROW_EXCEPTION
            ( 
               error::no_xml_tag() << error::line_number(l) 
            );
            result = false;
            break;
          }
          if(     ( _tags & tags::VALIDATE )
              and ( not bx::regex_match( match.prefix().first, match.prefix().second, empty ) ) )
          {
            size_t const l1( _line + std::count( ofirst, match.prefix().first, '\n' ) );
            size_t const l2( _line + std::count( ofirst, match.prefix().second, '\n' ) );
            if( l1 == l2 )
            { 
              BOOST_THROW_EXCEPTION
              ( 
                error::no_xml_tag() << error::line_number(l1)
                                    << error::section_name(match.prefix().str()) 
              );
            } 
            else 
            { 
              BOOST_THROW_EXCEPTION
              ( 
                error::no_xml_tag() << error::line_number(l1)
                                    << error::line_number(l2)
                                    << error::section_name(match.prefix().str()) 
              );
            } 
            result = false;
          }
          std::string const tag_name( match[1].str() );
           
          std::string::const_iterator const tag_end( std::find( match[0].second, _last, '>' ) );
          std::string::const_iterator const tag_begin( match[0].first );
          if( tag_end == _last )
          {
            if( not ( _tags & tags::VALIDATE ) ) break;
            size_t const l( _line + std::count( ofirst, tag_begin, '\n' ) );
            BOOST_THROW_EXCEPTION
            ( 
              error::no_end_tag() << error::line_number(l)
                                  << error::section_name(tag_name)
            );
            result = false;
            break;
          }
     
          // check options.
          tree::Base::t_Subsections::value_type subsection( tag_name );
          bool const no_subsections( *(tag_end-1) == '/' );
          result &= parse_options
                    ( 
                      subsection,
                      match[0].second,
                      no_subsections ? tag_end-1: tag_end, 
                      _line, _tags, tag_name
                    );
          
          if( no_subsections ) 
          {
            _first = tag_end + 1;
            if( _tags & tags::READ ) _base.push_back( subsection );
            continue;
          }
     
          bx::sregex const end_tag = bx::as_xpr("</") >> *bx::_s >> tag_name >> *bx::_s >> ">";
          if( not bx::regex_search( tag_end + 1, _last, match, end_tag ) )
          {
            if( not (_tags & tags::VALIDATE) ) continue;
            size_t const l( _line + std::count( ofirst, tag_begin, '\n' ) );
            BOOST_THROW_EXCEPTION
            ( 
              error::no_end_tag() << error::line_number(l)
                                  << error::section_name(tag_name)
            );
            result = false;
            _first = tag_end + 1;
            continue;
          }
          result &= parse( subsection, tag_end+1, match[0].first, _tags, 
                           _line+std::count(ofirst, tag_end+1, '\n') );
          _first = match[0].second;
          if( _tags & tags::READ ) _base.push_back( subsection );
        }
     
        return result;
      };

      boost::shared_ptr<tree::Base> parse( std::string const &_text )
      {
        boost::shared_ptr<tree::Base> result( new tree::Base );
        if( not parse( *result, _text.begin(), _text.end(), tags::READ_N_VALIDATE, 0 ) )
          return boost::shared_ptr<tree::Base>();
        return result;
      }
 
    } // namespace xml
  } // namespace load_n_save
} // namespace LaDa

