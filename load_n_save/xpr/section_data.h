#ifndef LADA_LNS_XPR_SECTIONDATA_H
#define LADA_LNS_XPR_SECTIONDATA_H

#include "LaDaConfig.h"

#include <boost/lexical_cast.hpp>
#include <boost/type_traits/remove_const.hpp>

#include "../tags.h"
#include "../access.h"
#include "../string_type.h"
#include "../parser_base.h"
#include "../exceptions.h"

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
        tags tag;
        //! Help of the section.
        t_String help;
      };

      namespace details
      {
        //! Virtual base class for section data.
        template<class T_SECTION>
          struct SectionDataBase
          {
            //! Virtual Destructor.
            virtual ~SectionDataBase() {}
            //! Returns a print-out string.
            virtual t_String print(t_String const &_depth, t_String const &_tab) const = 0;
            //! Parsing double dispatch.
            virtual bool parse( parser_base::Section const&,
                                T_SECTION const& _sec, version_type _version ) const = 0;
            //! Returns tag of this section.
            virtual tags tag() const { return unavailable; }
          };

        //! Holds section data.
        template< class T_SECTION, class T_DATA > class SectionData;

        //! Holds external type data.
        template<class T_SECTION, class T_DATA> 
          class SectionData< T_SECTION, boost::shared_ptr<T_DATA> > : public SectionDataBase<T_SECTION>
          {
            public:
        
              //! Virtual Destructor.
              virtual ~SectionData() {};
              //! Constructor.
              SectionData( boost::shared_ptr<T_DATA> _data) : data_(_data) {};
              //! Copy Constructor.
              SectionData   (const SectionData &_c)
                          : SectionDataBase<T_SECTION>(_c), data_(_c.data_) {};
              //! Returns a print-out string.
              virtual t_String print(t_String const &_depth, t_String const &_tab) const
                { return "External shared pointer.\n"; }
              //! Parsing double dispatch.
              virtual bool parse( parser_base::Section const& _parser, 
                                  T_SECTION const &, version_type _version ) const
              {
                if(not data_) BOOST_THROW_EXCEPTION(error::empty_shared_ptr());
                return LaDa::load_n_save::lns_access_adl( _parser, *data_, _version ); 
              }
            private:
              //! The data itself
              boost::shared_ptr<T_DATA> data_;
          };
      
        //! Holds external type data.
        template<class T_SECTION, class T_DATA> 
          class SectionData : public SectionDataBase<T_SECTION>
          {
            public:
      
              //! Virtual Destructor.
              virtual ~SectionData() {};
              //! Constructor.
              SectionData( T_DATA & _data) : data_(_data) {};
              //! Copy Constructor.
              SectionData  (const SectionData &_c)
                          : SectionDataBase<T_SECTION>(_c), data_(_c.data_) {};
              //! Returns a print-out string.
              virtual t_String print(t_String const &_depth, t_String const &_tab) const
                { return "External type.\n"; }
              //! Parsing double dispatch.
              virtual bool parse( parser_base::Section const& _parser, 
                                  T_SECTION const &, version_type _version ) const
                { return LaDa::load_n_save::lns_access_adl( _parser, data_, _version ); }
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
              //! Returns a print-out string.
              virtual t_String print(t_String const &_depth, t_String const &_tab) const
              {
                t_String tag = boost::lexical_cast<t_String>(data_.tag);
                return   _depth + "section " + data_.name + "\n"
                       + _depth + _tab + "help: " + data_.help + "\n"
                       + _depth + _tab + "tag:  " + tag + "\n";
              }
              //! Parsing double dispatch.
              virtual bool parse( parser_base::Section const& _parser,
                                  T_SECTION const &this_, version_type ) const
                { return _parser.regular(this_, data_); }
              //! Returns tag of this section.
              tags tag() const { return data_.tag; }
            private:
              //! The data itself
              regular_data data_;
          };

        } // namespace details.
    } // namespace xpr
  } // namespace load_n_save
} // namespace LaDa

#endif