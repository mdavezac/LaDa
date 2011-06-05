#ifndef LADA_LOADNSAVE_ACTION_SHAREDPTR_H
#define LADA_LOADNSAVE_ACTION_SHAREDPTR_H

#include "LaDaConfig.h"

#include <boost/shared_ptr.hpp>

#include "../xpr/section.h"
#include "../string_type.h"
#include "../access.h"
#include "action_base.h"

namespace LaDa 
{
  namespace load_n_save
  {
    template<class T_ARCHIVE, class T_TYPE>
      bool lns_access(T_ARCHIVE &_ar, boost::shared_ptr<T_TYPE> &_ptr, version_type const _version)
      {
        xpr::Section result;
        result.set_data(*_ptr);
        return _ar & result;
      };

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
