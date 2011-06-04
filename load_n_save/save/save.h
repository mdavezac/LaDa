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
    namespace save
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
          boost::shared_ptr<tree::Base> operator()(xpr::Section const& _sec ) const;
     
          //! Loads an archive.
          bool is_loading() const { return false; } 

        protected:
          //! Class for parsing a section.
          class Section;
      };

      class Save :: Section : public parser_base::Section
      {
        public:
          //! Constructor.
          Section(boost::shared_ptr<tree::Section> const &_tree) : tree_(_tree) {}
          //! Copy constructor.
          Section(Section & _c) : tree_(_c.tree_) {}

          bool operator()(xpr::Section const& _sec) const
            { return operator&(_sec); }
          virtual bool operator&( xpr::Section const& _sec ) const
            { return _sec.parse(*this); }

          //! Saves an archive.
          bool is_loading() const { return false; } 

          //! Parses content.
          virtual bool content( xpr::Section const & _sec, t_String const &_name = "" ) const;
          //! Double dispatch.
          virtual bool regular( xpr::Section const &_sec, xpr::regular_data const&_data ) const;
          //! Returns the current tree;
          boost::shared_ptr<tree::Section> tree() const { return tree_; }

        protected:
          //! Parses content of an expression range;
          template<class T> bool content_(T const &_range) const;
          //! Reference to the tree being built.
          boost::shared_ptr<tree::Section> tree_;
      };

    } // namespace initializer.
  } // namespace load_n_save
} // namespace LaDa


#endif
