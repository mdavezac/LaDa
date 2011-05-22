#include "LaDaConfig.h"

#include <boost/iterator/filter_iterator.hpp>

#include "save.h"
#include "../sequencer/binary_range.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace save
    {
      boost::shared_ptr<tree::Base> Save::operator()(xpr::Section const& _sec ) const
      { 
        boost::shared_ptr<tree::Base> result(new tree::Base);
        boost::shared_ptr<tree::Section> text(new tree::Section(""));
        Save::Section section(text); section(_sec);
        result->push_back(*text);
        return result;
      }

      // by defining these functions this way, we won't have section.h depending on sequencer_range.h.
      // makes for faster compilation. These functions are needed for
      // recurrence over options linked by "and" and ""or" operators.
      typedef sequencer::BinaryRange
              <
                xpr::Section::const_iterator::subsection,
                xpr::Section::const_iterator::option
              > OpAndRange;
      template<class T>
        bool Save::Section :: content_(T const &_range) const
        {
          bool error = false;
          OpAndRange::value_type i_first(_range.first);
          OpAndRange::value_type const &i_end(_range.second);
          for(; i_first != i_end; ++i_first )
            if( i_first.is_first() )
            {
              boost::shared_ptr<tree::Section> treesec(new tree::Section(""));
              Save::Section section(treesec);
              if( not i_first.first()->parse(section) ) return false;
              treesec->copy_to(*tree_);
            }
            else if( i_first.is_second() )
              tree_->push_back( i_first.second()->name, i_first.second()->str() );
            else if( i_first.is_start_group() )
            { 
              OpAndRange::value_type i_group(_range.first);
              find_group_end( i_group, i_end );
              if( not content_(OpAndRange(++i_first, i_group)) ) error = true;
              i_first = i_group;
            }
          return error;
        }

      bool Save::Section :: regular( xpr::Section const & _sec, xpr::regular_data const &_data ) const
      {
        tree_->name = _data.name;
        tree_->fresh = true;

        OpAndRange const range
        (
          _sec.sequence().begin(), _sec.sequence().end(),
          _sec.subsections_begin(), _sec.subsections_end(),
          _sec.options_begin(), _sec.options_end() 
        );
        return content_(range);
      }


    } // namespace save.
  } // namespace load_n_save
} // namespace LaDa


