#ifndef LADA_LOADNSAVE_SAVE_SAVE_H
#define LADA_LOADNSAVE_SAVE_SAVE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../tree/tree.h"
#include "../xpr/section.h"
#include "../string_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace load
    {
      //! Links text with operators for loading purposes.
      class Save 
      {
        public:
          //! Verbose or quiet output.
          bool verbose;
     
          //! Constructor.
          Save() : verbose(true) {}
          //! Copy Constructor.
          Save( Save const &_c ) : verbose(_c.verbose) {}

          //! Returns an in put/output tree.
          boost::shared_ptr<tree::Section> operator()(xpr::Section const& _sec ) const
            { return Section()(_sec).tree(); }
     
          //! Loads an archive.
          bool is_loading() const { return false; } 

        protected:
          //! Class for parsing a section.
          class Section;
          //! Class for parsing an option.
          class Option;
      };

      class Load :: Section : public parser_base::Section
      {
        public:
          //! Constructor.
          Section() : tree_(new tree::Section) {}
          //! Copy constructor.
          Section(Section & _c) : tree_(_c.tree) {}

          bool operator()(xpr::Section const& _sec) const
            { return operator&(_sec); }
          virtual bool operator&( xpr::Section const& _sec ) const
            { return _sec.parse(*this); }

          //! Loads an archive.
          bool is_loading() const { return false; } 

          //! Parses content.
          virtual bool content( xpr::Section const & _sec, t_String const &_name = "" ) const;
          //! Double dispatch.
          virtual bool regular( xpr::Section const &_sec, xpr::regular_data const&_data ) const;
          //! Returns the current tree;
          boost::shared_ptr<tree::Section> tree() const { return tree_; }

        protected:
          //! Reference to the tree being built.
          boost::shared_ptr<tree::Section> tree_;
      };

      class Load :: Option 
      {
          friend class Section;
        public:
          //! Parsing Constructor.
          Option( Section const &_c) : tree_(_c.tree_) {}
          //! Copy Constructor.
          Option( Option const &_c ) : tree_(_c.tree_) {}

          //! Starts parsing.
          bool operator()( t_String const& _name, xpr::Option const& _sec ) const;

        protected:
          //! Checks a value against a regex.
          bool parse_( t_String const& _regex, t_String const &_value ) const;
          //! Xml-like tree.
          boost::shared_ptr<tree::Section> tree_;
      };

    } // namespace initializer.
  } // namespace load_n_save
} // namespace LaDa


#endif
