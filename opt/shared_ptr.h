#ifndef LADA_OPT_SHARED_PTR_H
#define LADA_OPT_SHARED_PTR_H

#include "LaDaConfig.h"

#include <boost/type_traits/is_array.hpp>

namespace LaDa 
{

  template<class T> 
    class SharedPtr
    {
      public:
        //! Null constructor.
        SharedPtr() : previous_(NULL), next_(NULL), object_(NULL) {}
        //! Ownership constructor.
        SharedPtr(T * const _object) : previous_(NULL), next_(NULL), object_(_object) {}
        //! Shared ownership constructor.
        SharedPtr(SharedPtr<T> &_t) : previous_(NULL), next_(NULL), object_(NULL)
        {
          //! other object is NULL
          if(_t.object_ == NULL) return;
          previous_ = &_t;
          next_ = _t.next_;
          _t.next_ = this;
          object = _t.object;
        }
        //! Destruction
        virtual ~SharedPtr()
        {
          if(next_ != NULL) next_->previous_ = previous_;
          if(previous_ != NULL) previous_->next_ = next_;
          if(object_ != NULL and next_ == NULL and previous_ == NULL)
            deallocate_object_(boost::is_array<T>());
          previous_ = NULL;
          next_ = NULL;
          object_ = NULL;
        }




      private:
        void deallocate_object(T *


        //! Pointer to previous owner in list of owners.
        SharedPtr<T> *previous_;
        //! Pointer to next owner in list of owners.
        SharedPtr<T> *next_;
        //! Pointee.
        T * object_;
    };

} // namespace LaDa

#endif
