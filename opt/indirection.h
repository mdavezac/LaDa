//
//  Version: $Id$
//
#ifndef _OPT_INDIRECTION_H_
#define _OPT_INDIRECTION_H_

#include<boost/shared_ptr.hpp>

namespace Indirection 
{
  //! Base class for wrapping objects.
  template< class T_OBJECT >
    class Object
    {
      public:
        //! Type of the matrix.
        typedef T_OBJECT t_Object;
  
        //! Constructor.
        Object() {}
        //! Copy Constructor.
        Object( const Object &_c ) : object_(_c.object_) {} 
        //! Returns a reference to the object.
        operator t_Object& () { return object_; }
        //! Returns a reference to the object.
        operator const t_Object& () const { return object_; }
        //! Returns a reference to the object.
        t_Object& operator()() { return object_; }
        //! Returns a reference to the object.
        const t_Object& operator()() const { return object_; }
  
      protected:
        //! The coefficients themselves.
        t_Object object_;
    };
  //! Base class for wrapping pointer and interface with them as objects.
  template< class T_OBJECT >
    class Pointer
    {
      public:
        //! Type of the matrix.
        typedef T_OBJECT t_Object;
   
        //! Constructor.
        Pointer() : object_( new t_Object ) {}
        //! Copy Constructor.
        Pointer( const Pointer &_c ) : object_(_c.object_) {} 
        //! Returns a reference to the object.
        operator t_Object& () { return *object_; }
        //! Returns a reference to the object.
        operator const t_Object& () const { return *object_; }
        //! Returns a reference to the object.
        t_Object& operator()() { return *object_; }
        //! Returns a reference to the object.
        const t_Object& operator()() const { return *object_; }
   
      protected:
        //! The coefficients themselves.
        boost::shared_ptr<t_Object> object_;
    };
}
#endif
