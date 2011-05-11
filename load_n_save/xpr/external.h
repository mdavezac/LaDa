//
//  Version: $Id: external.h 1250 2009-07-26 21:04:07Z davezac $
//

#ifndef _LADA_LNS_XPR_SECTIONDATA_H_
#define _LADA_LNS_XPR_SECTIONDATA_H_

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/lexical_cast.hpp>

#include "../string_type.h"
#include "../parser_base/section.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      //! Regular section data.
      struct regular_data 
      {
        //! Name of the section.
        t_String name;
        //! Tags of the section.
        size_t tag;
        //! Help of the section.
        t_String help;
      };

      namespace details
      {
        namespace tracking = boost::xpressive::detail;

        //! Virtual base class for section data.
        template<class T_SECTION>
          struct SectionDataBase
          {
            //! Virtual Destructor.
            virtual ~SectionDataBase() {}
            //! Returns name.
            virtual t_String const& name() const = 0;
            //! Returns tag.
            virtual size_t tag() const = 0;
            //! Returns help.
            virtual t_String const& help() const = 0;
            //! Returns a print-out string.
            virtual t_String print(t_String const &_depth, t_String const &_tab) const = 0;
            //! Parsing double dispatch.
            virtual bool parse( parser_base::Section const&, T_SECTION const& _sec ) const = 0;
          };

        //! Holds section data.
        template< class T_SECTION, class T_DATA = regular_data > class SectionData;
      
        //! Holds external type data.
        template<class T_SECTION, class T_DATA> 
          class SectionData<T_DATA, T_DATA> : public SectionDataBase<T_SECTION>
          {
            public:
      
              //! Virtual Destructor.
              virtual ~ExternalData() {};
              //! Constructor.
              ExternalData( regular_data const& _data) : data_(_data) {};
              //! Copy Constructor.
              ExternalData   (const ExternalData &_c)
                           : SectionDataBase<T_SECTION>(_c), data_(_c.data_) {};
              //! Returns name.
              virtual t_String const& name() const { return "Should not be here."; }
              //! Returns tag.
              virtual size_t tag() const { return 0; }
              //! Returns help.
              virtual t_String const& help() const { return "Should not be here."; }
              //! Returns a print-out string.
              virtual t_String print(t_String const &_depth, t_String const &_tab) const
                { return "External type.\n" }
              //! Parsing double dispatch.
              virtual bool parse( parser_base::Section const& _parser, T_SECTION const & ) const
                { return acess<T_DATA>()( _parser, data_ ); }
            private:
              //! The data itself
              T_DATA &data_;
          };

        //! Holds regular data.
        template<class T_SECTION> 
          class SectionData<T_SECTION, regular_data> : public SectionDataBase<T_SECTION>
          {
            public:
  
              //! Virtual Destructor.
              virtual ~SectionData() {};
              //! Constructor.
              SectionData() { data_.name = ""; data_.tag = 0; data_.help = ""; }
              //! Constructor.
              SectionData( regular_data const& _data) : data_(_data) {};
              //! Copy Constructor.
              SectionData   (const SectionData &_c)
                          : SectionDataBase<T_SECTION>(_c), data_(_c.data_) {};
              //! Returns name.
              virtual t_String const& name() const { return data_.name; }
              //! Returns tag.
              virtual size_t tag() const { return data_.tag; }
              //! Returns help.
              virtual t_String const& help() const { return data_.help; }
              //! Returns a print-out string.
              virtual t_String print(t_String const &_depth, t_String const &_tab) const
              {
                t_String tag = boost::lexical_cast<t_String>(data_.tag);
                return   _depth + "section " + data_.name + "\n"
                       + _depth + _tab + "help: " + data_.help + "\n"
                       + _depth + _tab + "tag:  " + tag + "\n";
              }
              //! Parsing double dispatch.
              virtual bool parse( parser_base::Section const& _parser, T_SECTION const &this_ ) const
                { return _parser.regular(this_, data_); }
            private:
              //! The data itself
              regular_data data_;
          };

        } // namespace details.
    } // namespace xpr
  } // namespace load_n_save
} // namespace LaDa

#endif
