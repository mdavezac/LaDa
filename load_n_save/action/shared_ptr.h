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
    template<class T_TYPE>
      struct lns_access< boost::shared_ptr<T_TYPE> > 
      {
        template<class T_ARCHIVE>
          bool operator()(T_ARCHIVE const &_ar, boost::shared_ptr<T_TYPE> const &_var)
          {
            xpr::Section result;
            result.set_data(*_var);
            return _ar & result;
          }
      };
//   namespace action_
//   {
//     //! Wrapper around shared ptr.
//     template< class T_TYPE >
//       class SharedPtr
//       {
//         public:
//           //! The map type
//           typedef boost:shared_ptr<T_TYPE> t_ptr;
//           //! \brief Constructor. 
//           //! \params _var
//           SharedPtr (t_ptr const &_var) : var_(_var) {}
//           //! CopyConstructor. 
//           SharedPtr(SharedPtr const &_c) : var_(_c.var_) {}
//
//           template<class T_ARCHIVE>
//             bool lns_access(T_ARCHIVE const &_ar) { return _ar & *var_; }
//         protected:
//           //! Holds reference to variable.
//           t_ptr var_;
//       };
//    
//   } // namespace action.

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
