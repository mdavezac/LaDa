#include "LaDaConfig.h"

#include <boost/iterator/filter_iterator.hpp>

#include "load.h"
#include "../tags.h"
#include "../sequencer/binary_range.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace load
    {

      namespace details
      {
        // by defining these functions this way, we won't have section.h depending on sequencer_range.h.
        // makes for faster compilation. These functions are needed for
        // recurrence over options linked by "and" and ""or" operators.
        typedef sequencer::BinaryRange
                <
                  xpr::Section::const_iterator::subsection,
                  xpr::Section::const_iterator::option
                > OpAndRange;
        struct content_impl
        {
          typedef std::pair<bool, bool> t_ReturnId;
          t_ReturnId static id( Load::Section const& _self,
                                OpAndRange const&_range,
                                t_String const &_name );
          t_ReturnId static id2( Load::Section const& _self, 
                                 OpAndRange const &_range,
                                 t_String const &_name);
          bool static content( Load::Section const& _self, OpAndRange const&_range, 
                               t_String const &_name, bool _grammar_only);
          bool static content2( Load::Section const& _self, OpAndRange const &_range, 
                                t_String const &_name, bool _grammar_only);
        };
        // Allows iteration over fresh (previously unparsed) options and sections.
        struct Fresh
        {
          bool operator()( tree::Section const& _sec ) const { return _sec.fresh; }
          bool operator()( tree::Option const& _sec ) const { return _sec.fresh; }
        };
        
        // Prints out fresh options and sequence on the assumption they were not parsed.
        void print_unknowns( std::string const& _name, tree::Section const& _sec );
      }


      bool Load::Section :: content( xpr::Section const & _sec, t_String const &_name ) const
      {
        details::OpAndRange range
        (
          _sec.sequence().begin(), _sec.sequence().end(),
          _sec.subsections_begin(), _sec.subsections_end(),
          _sec.options_begin(), _sec.options_end() 
        );

        bool found = details::content_impl::content( *this, range, _name, grammar_only_ );
        if( (not found) and verbose_ and (not grammar_only_) )
        {
          const std::string dummy( _name.size() > 0 ? " of section " + _name: " " );
          std::cerr << "Some subsections" + dummy + " were not found "
                       "or did not parse.\n";
        }
        return found;
      }

      bool Load::Section :: regular( xpr::Section const & _sec, xpr::regular_data const &_data ) const
      {
        tree::Section::const_iterator :: subsection i_first( tree_.subsections_begin(_data.name) );
        tree::Section::const_iterator :: subsection const i_last( tree_.subsections_end(_data.name) );

        bool found(false);
        bool parse_error(false);
        for(; i_first != i_last; ++i_first )
        {
          if( not i_first->fresh ) continue;

          Section parser(*this, *i_first, version_);
          if( not parser.id_options_(_data.name, _sec) ) continue;
          
          found = true;
          if( grammar_only_ ) break;

          // parses section.
          parser.grammar_only_ = false;
          if( not parser.content(_sec, _data.name) ) parse_error = true;
          i_first->fresh = false; // marks section as read.
          details::print_unknowns( _data.name, *i_first );

          break;
        }
        

        if( found )
        {
          if( verbose_ and parse_error )
            std::cerr << "Found section " + _data.name + " but could not parse it.\n";
          return not parse_error;
        }
        if( _data.tag & load_n_save::required)
        {
          if( verbose_ and (not grammar_only_) )
            std::cerr << "Did not find required section " + _data.name + ".\n";
          return false;
        }

        return true;
      }

      bool Load::Section :: id_options_(t_String const &_name, xpr::Section const& _sec ) const
      {
        details::OpAndRange range
        (
          _sec.sequence().begin(), _sec.sequence().end(),
          _sec.subsections_begin(), _sec.subsections_end(),
          _sec.options_begin(), _sec.options_end() 
        );

        details::content_impl::t_ReturnId found = details::content_impl::id( *this, range, _name );
        if( found.first or (not found.second) ) return true;
        return false;
      }



      namespace details
      {
        bool content_impl::content( Load::Section const& _self, OpAndRange const &_range, 
                                    t_String const &_name, bool _grammar_only)
        {
          bool found = false;
          OpAndRange range(_range);
          for(; (bool)range; ++range )
          {
            found = content2(_self, range, _name, true);
            if( found ) break;
          }
          // not found, or not parsing now.
          if( (not found) or _grammar_only ) return found;
        
          // parses.
          return content2(_self, range, _name, _grammar_only);
        };
        bool content_impl::content2( Load::Section const& _self, OpAndRange const &_range,
                                     t_String const &_name, bool _grammar_only) 
        {
          OpAndRange::value_type i_first(_range.first);
          OpAndRange::value_type const &i_end(_range.second);
          bool error = false;
          for(; i_first != i_end; ++i_first )
            if( i_first.is_first() )
            {
              Load::Section section(_self);
              section.grammar_only_ = _grammar_only;
              if( not i_first.first()->parse(section, section.version_) ) return false;
            }
            else if( i_first.is_second() )
            {
              Load::Option option(_self);
              option.grammar_only_ = _grammar_only;
              if( option(_name, *i_first.second() ) ) continue;
              // if here, "error" occured. 
              if(     i_first.second()->tag & load_n_save::required
                  and  _self.verbose_
                  and (not _grammar_only) )
                std::cerr << "Could not find required option " << i_first.second()->name 
                          << " in section " + _name + ".\n";
              return false;
            }
            else if( i_first.is_start_group() )
            { 
              OpAndRange::value_type i_group(_range.first);
              find_group_end( i_group, i_end );
              bool const found = content( _self, OpAndRange(++i_first, i_group), _name, _grammar_only );
              if( not found ) error = true;
              i_first = i_group;
            }
#           ifdef _LADADEBUG
              else { __THROW_ERROR("Should not be here.\n"); }
#           endif
          return  not error;
        }

        content_impl::t_ReturnId content_impl::id( Load::Section const& _self,
                                                   OpAndRange const &_range,
                                                   t_String const &_name )
        {
          t_ReturnId found(false, false);
          bool possible_noid = false; // is true if one sequence of options has no id.
          OpAndRange range(_range);
          for(; (bool)range; ++range )
          {
            found = id2(_self, range, _name);
            if( found.first and found.second ) return found;
            if( not found.second ) possible_noid = true;
          }
          if( possible_noid ) found.second = false;
          return found;
        };

        content_impl::t_ReturnId content_impl::id2( Load::Section const& _self,
                                                    OpAndRange const &_range, 
                                                    t_String const &_name) 
        {
          OpAndRange::value_type i_first(_range.first);
          OpAndRange::value_type const &i_end(_range.second);
          bool has_id(false);
          for(; i_first != i_end; ++i_first )
            if( i_first.is_second() )
            {
              if( not (i_first.second()->tag & load_n_save::id) ) continue;
              has_id = true;
              Load::Option option(_self); option.grammar_only_ = true;
              bool const found = option(_name, *i_first.second() );
              if( not found ) return t_ReturnId(false, true);
            }
            else if( i_first.is_start_group() )
            { 
              OpAndRange::value_type i_group(_range.first);
              find_group_end( i_group, i_end );
              t_ReturnId const found = id( _self, OpAndRange(++i_first, i_group), _name );

              if( (not found.first) and found.second ) return t_ReturnId(false, true);

              i_first = i_group;
              if( found.second ) has_id = true;
            }
#           ifdef _LADADEBUG
              else if( i_first.is_first() ) continue;
              else { __THROW_ERROR("Should not be here.\n"); }
#           endif
          return  t_ReturnId(has_id,has_id);
        }

        // Prints out fresh options and sequence on the assumption they were not parsed.
        void print_unknowns( std::string const& _name, tree::Section const& _sec )
        {
          typedef boost::filter_iterator<Fresh, tree::Section::const_iterator::subsection > t_itersec;
          t_itersec i_sec( _sec.subsections_begin(), _sec.subsections_end() );
          t_itersec const i_sec_end( Fresh(), _sec.subsections_end(), _sec.subsections_end() );
          if( i_sec != i_sec_end )
          {
            std::cerr << "Unknown sections in " + _name + ":";
            for(; i_sec != i_sec_end; ++i_sec )
              std::cerr << " " << i_sec->name;
            std::cerr << "\n";
          }
        
          typedef boost::filter_iterator<Fresh, tree::Section::const_iterator::option > t_iterop;
          t_iterop i_op( _sec.options_begin(), _sec.options_end() );
          t_iterop const i_op_end( _sec.options_end(), _sec.options_end() );
          if( i_op != i_op_end )
          {
            std::cerr << "Unknown options in " + _name + ":";
            for(; i_op != i_op_end; ++i_op )
              std::cerr << " " << i_op->name;
            std::cerr << "\n";
          }
        }
      }

    } // namespace load.
  } // namespace load_n_save
} // namespace LaDa

