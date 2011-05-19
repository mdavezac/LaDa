#ifndef LADA_LOADNSAVE_LOAD_LOAD_H
#define LADA_LOADNSAVE_LOAD_LOAD_H

#include "LaDaConfig.h"

#include "../tree/tree.h"
#include "../xpr/section.h"
#include "../string_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace load
    {
      namespace details
      {
        struct content_impl;
      }

      //! Links text with operators for loading purposes.
      class Load 
      {
          friend class details::content_impl;
          //! Class parsing sections.
          class Section;
          //! Class parsing options.
          class Option;
        public:
          //! Verbose or quiet output.
          bool verbose;
     
          //! Constructor.
          Load() : verbose(true) {}
          //! Copy Constructor.
          Load( Load const &_c ) : verbose(_c.verbose) {}

          //! Starts parsing.
          bool operator()( tree::Base const& _tree, xpr::Section const& _sec ) const;
          //! Starts parsing.
          bool operator()( tree::Section const& _tree, xpr::Section const& _sec ) const;
     
          //! Loads an archive.
          bool is_loading() const { return true; } 
      };

      class Load :: Section : public parser_base::Section
      {
          friend class details::content_impl;
          friend class Option;
        public:
          //! Constructor.
          Section   ( tree::Section const& _tree, bool _verbose ) 
                  : tree_(_tree), verbose_(_verbose), grammar_only_(false) {}
          //! Parsing Constructor.
          Section   ( Section const &_c, tree::Section const& _tree ) 
                  : tree_(_tree), verbose_(_c.verbose_), grammar_only_(true) {}
          //! Copy Constructor.
          Section   ( Section const &_c )
                  : tree_(_c.tree_), verbose_(_c.verbose_), grammar_only_(_c.grammar_only_) {}

          //! Starts parsing.
          bool operator()( xpr::Section const& _sec ) const
            { return operator&(_sec); }

          //! Starts parsing.
          virtual bool operator&( xpr::Section const& _sec ) const
            { return _sec.parse(*this); }

          //! Loads an archive.
          bool is_loading() const { return true; } 

          //! Parses content.
          virtual bool content( xpr::Section const & _sec, t_String const &_name = "" ) const;
          //! Double dispatch.
          virtual bool regular( xpr::Section const &_sec, xpr::regular_data const&_data ) const;

        protected:
          //! Checks id_options.
          bool id_options_(t_String const& _name, xpr::Section const& _sec ) const;

          //! Xml-like tree.
          tree::Section const &tree_;
          //! Verbose or quiet output.
          bool verbose_;
          //! Whether to parse options or simply check syntax.
          bool grammar_only_;
      };
     
      class Load :: Option 
      {
          friend class Section;
          friend class details::content_impl;
        public:
          //! Parsing Constructor.
          Option   ( Section const &_c )
                 : tree_(_c.tree_), verbose_(_c.verbose_), grammar_only_(_c.grammar_only_) {}
          //! Copy Constructor.
          Option   ( Option const &_c )
                 : tree_(_c.tree_), verbose_(_c.verbose_), grammar_only_(_c.grammar_only_) {}

          //! Starts parsing.
          bool operator()( t_String const& _name, xpr::Option const& _sec ) const;

        protected:
          //! Checks a value against a regex.
          bool parse_( t_String const& _regex, t_String const &_value ) const;
          //! Xml-like tree.
          tree::Section const &tree_;
          //! Verbose or quiet output.
          bool verbose_;
          //! Whether to parse options or simply check syntax.
          bool grammar_only_;
      };

      inline bool Load::operator()( tree::Base const& _tree, xpr::Section const& _sec ) const
      {
        tree::Base base(_tree);
        tree::Section text("");
        text.tree::Base::swap(base);
        return Section(text, verbose)(_sec); 
      }
      inline bool Load::operator()( tree::Section const& _tree, xpr::Section const& _sec ) const
        { return Section(_tree, verbose)(_sec); }
    } // namespace initializer.
  } // namespace load_n_save
} // namespace LaDa


#endif
