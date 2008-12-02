//
//  Version: $Id: base.h 844 2008-11-08 01:22:54Z davezac $
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <opt/types.h>
#include <opt/debug.h>

namespace LaDa
{
  namespace Print
  {
    void format_description(std::ostream& _stream,
                            const std::string& desc, 
                            unsigned first_column_width,
                            unsigned line_length);
    //! \brief Given a string 'par', that contains no newline characters
    //!        outputs it to '_stream' with wordwrapping, that is, as several line.
    //! \detail Each output line starts with 'indent' space characters,
    //!         following by characters from 'par'. The total length of line is
    //!         no longer than 'line_length'.
    void format_paragraph(std::ostream& _stream,
                          std::string par,
                          unsigned indent,
                          unsigned line_length);

    void format_one( std::ostream& _stream, const std::string& _first, 
                     const std::string& _second,
                     size_t _first_column_width, size_t _line_length)
    {
      // Don't use since g++ 2.96 is buggy on it.
      _stream << _first;

      if (!_second.empty())
      {
        for(unsigned pad = _first_column_width - _first.size();  
            pad > 0; --pad) _stream.put(' ');
      
        format_description(_stream, _second, _first_column_width, _line_length);
      }
    }

    void format_description(std::ostream& _stream,
                            const std::string& desc, 
                            unsigned first_column_width,
                            unsigned line_length)
    {
      // we need to use one char less per line to work correctly if actual
      // console has longer lines
      assert(line_length > 1);
      if (line_length > 1)
      {
          --line_length;
      }

      // line_length must be larger than first_column_width
      // this assert may fail due to user error or environment conditions!
      assert(line_length > first_column_width);

      // Note: can't use 'tokenizer' as name of typedef -- borland
      // will consider uses of 'tokenizer' below as uses of
      // boost::tokenizer, not typedef.
      typedef boost::tokenizer<boost::char_separator<char> > tok;
      
      tok paragraphs( desc, boost::char_separator<char>("\n", "", boost::keep_empty_tokens));
      
      tok::const_iterator       par_iter = paragraphs.begin();                
      const tok::const_iterator par_end = paragraphs.end();

      while (par_iter != par_end)  // paragraphs
      {
        format_paragraph(_stream, *par_iter, first_column_width, 
                         line_length);
        
        ++par_iter;
        
        // prepair next line if any
        if (par_iter != par_end)
        {
          _stream << '\n';
        
          for(unsigned pad = first_column_width; pad > 0; --pad)
            _stream.put(' ');
        }            
      }  // paragraphs
    }
    
    void format_paragraph(std::ostream& _stream,
                          std::string par,
                          unsigned indent,
                          unsigned line_length)
    {                    
      // Through reminder of this function, 'line_length' will
      // be the length available for characters, not including
      // indent.
      assert(indent < line_length);
      line_length -= indent;

      // index of tab (if present) is used as additional indent relative
      // to first_column_width if paragrapth is spanned over multiple
      // lines if tab is not on first line it is ignored
      std::string::size_type par_indent = par.find('\t');

      if (par_indent == std::string::npos) par_indent = 0;
      else
      {
        // only one tab per paragraph allowed
        __ASSERT(count(par.begin(), par.end(), '\t') > 1,
                 "Only one tab per paragraph is allowed")
      
        // erase tab from string
        par.erase(par_indent, 1);

        // this assert may fail due to user error or 
        // environment conditions!
        assert(par_indent < line_length);

        // ignore tab if not on first line
        if (par_indent >= line_length) par_indent = 0;
      }
      
      if (par.size() < line_length) _stream << par;
      else
      {
        std::string::const_iterator       line_begin = par.begin();
        const std::string::const_iterator par_end = par.end();

        bool first_line = true; // of current paragraph!        
      
        while (line_begin < par_end)  // paragraph lines
        {
          if (!first_line)
          {
            // If line starts with space, but second character
            // is not space, remove the leading space.
            // We don't remove double spaces because those
            // might be intentianal.
            if ((*line_begin == ' ') &&
                ((line_begin + 1 < par_end) &&
                 (*(line_begin + 1) != ' ')))
            {
                line_begin += 1;  // line_begin != line_end
            }
          }

          // Take care to never increment the iterator past
          // the end, since MSVC 8.0 (brokenly), assumes that
          // doing that, even if no access happens, is a bug.
          unsigned remaining = distance(line_begin, par_end);
          std::string::const_iterator line_end = line_begin + 
              ((remaining < line_length) ? remaining : line_length);
      
          // prevent chopped words
          // Is line_end between two non-space characters?
          if ((*(line_end - 1) != ' ') &&
              ((line_end < par_end) && (*line_end != ' ')))
          {
            // find last ' ' in the second half of the current paragraph line
            std::string::const_iterator last_space =
                find(std::reverse_iterator<std::string::const_iterator>(line_end),
                     std::reverse_iterator<std::string::const_iterator>(line_begin),
                     ' ')
                .base();
        
            if (last_space != line_begin)
            {                 
              // is last_space within the second half ot the 
              // current line
              if ((unsigned)distance(last_space, line_end) < 
                  (line_length - indent) / 2)
              {
                  line_end = last_space;
              }
            }                                                
          } // prevent chopped words
       
          // write line to stream
          copy(line_begin, line_end, std::ostream_iterator<char>(_stream));
        
          if (first_line)
          {
            indent += par_indent;
            first_line = false;
          }

          // more lines to follow?
          if (line_end != par_end)
          {
            _stream << '\n';
        
            for(unsigned pad = indent; pad > 0; --pad)
              _stream.put(' ');
          }
        
          // next line starts after of this line
          line_begin = line_end;              
        } // paragraph lines
      }          
    }                              
        
  } // namespace Print
} // namespace LaDa

