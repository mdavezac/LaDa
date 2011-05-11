//
//  Version: $Id: section.cc 1266 2009-08-10 05:01:26Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/debug.h>
#include "section.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      void merge_options( Section &_a, Section const &_b )
      {
        for
        (
          Section::const_iterator::option
            i_first = _b.options_begin(), 
            i_end = _b.options_end();
          i_first != i_end;
          ++i_first 
        ) _a.push_back( *i_first );
      }
      void merge_subsections( Section &_a, Section const &_b )
      {
        for
        (
          Section::const_iterator::subsection
            i_first = _b.subsections_begin(), 
            i_end = _b.subsections_end();
          i_first != i_end;
          ++i_first 
        ) _a.push_back( *i_first );
      }
      void merge_all( Section &_a, Section const &_b )
      {
        merge_options( _a, _b );
        merge_subsections( _a, _b );
      }
      Section sum( Section _a, Section _b, bool _andop )
      {
        bool const a_incomplete( _a.incomplete() );
        bool const b_incomplete( _b.incomplete() );
        if( not ( a_incomplete and b_incomplete ) )
        {
          Section result;
          result.push_back(_a);
          result.push_back(_b);
          result.sequence() = _andop ? "11": "1|1";
          return result; 
        }
        else if( a_incomplete and b_incomplete ) 
        {
          merge_all( _a, _b );
          if( _andop )  _a.sequence() &= _b.sequence();
          else  _a.sequence() |= _b.sequence();
          return _a;
        }
        else
        {
          Section &c( a_incomplete ? _a: _b );
          c.push_back( a_incomplete ? _b: _a  );
          if( _andop ) c.sequence() &= "1";
          else         c.sequence() &= "1";
          return c;
        }
      };
      Section operator&&( Section _a, Section _b )
        { return sum(_a, _b, true); }
      Section operator||( Section _a, Section _b )
        { return sum(_a, _b, false); }
      
      Section sum( Section& _a, Option const& _b, bool _andop )      
      {
        if( _a.incomplete() )
        {
          _a.push_back( _b ); 
          if( _andop ) _a.sequence() &= "2";
          else         _a.sequence() |= "2";
          return _a;
        }
        else 
        {
          Section c;
          c.push_back(_a);
          c.push_back(_b);
          c.sequence() = _andop ? "12": "1|2";
          return c;
        }
      }


      Section operator&&( Section _a, Option const &_b ) { return sum(_a, _b, true); }
      Section operator&&( Option const &_b, Section _a ) { return sum(_a, _b, true); }
      Section operator||( Section _a, Option const &_b ) { return sum(_a, _b, false); }
      Section operator||( Option const &_b, Section _a ) { return sum(_a, _b, false); }


      Section sum( Option const& _a, Option const& _b, bool _andop )   
      {
        Section result;
        result.push_back(_a);
        result.push_back(_b);
        result.sequence() = _andop ? "22": "2|2";
        return result;
      }
      Section operator&&( Option const& _a, Option const &_b ) { return sum( _a, _b, true ); }
      Section operator||( Option const& _a, Option const &_b ) { return sum( _a, _b, false ); }
      
      //  Insert content.
      Section operator<<( Section _a, Section const& _b )
      {
        LADA_DOASSERT( not _a.incomplete(), "Inserting into incomplete section.\n" )
        merge_all( _a, _b );
        _a.sequence() &= _b.sequence();
        return _a;
      }

      //! Insert option.
      Section operator<<( Section _a, Option const& _b )
      {
        LADA_DOASSERT( not _a.incomplete(), "Inserting into incomplete section.\n" )
        _a.push_back(_b); 
        _a.sequence() &= "2";
        return _a;
      }


    } // namespace xpr.
  } // namespace load_n_save

} // namespace LaDa


