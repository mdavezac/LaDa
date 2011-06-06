#ifndef LADA_LOADNSAVE_XPR_SHAREDPTR_H
#define LADA_LOADNSAVE_XPR_SHAREDPTR_H

#include "LaDaConfig.h"

#include <boost/shared_ptr.hpp>

#include "../string_type.h"
#include "section_data.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      namespace details
      {
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
                { return LaDa::load_n_save::lns_access_adl( _parser, *data_, _version ); }
            private:
              //! The data itself
              boost::shared_ptr<T_DATA> data_;
          };
      }
    }


    //! Returns an Enum action.
    template<class T_TYPE >
      boost::shared_ptr<T_TYPE> copy(T_TYPE const &_a)
        { return boost::shared_ptr<T_TYPE>( new T_TYPE(_a) ); }
    //! Returns an Enum action.
    template<class T_TYPE >
      boost::shared_ptr<T_TYPE> copy(T_TYPE &_a)
        { return boost::shared_ptr<T_TYPE>( new T_TYPE(_a) ); }
  }
}

#endif
