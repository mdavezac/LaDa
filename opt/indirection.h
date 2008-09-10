//
//  Version: $Id$
//
#ifndef _OPT_INDIRECTION_H_
#define _OPT_INDIRECTION_H_

namespace opt 
{
  //! Pure base class for wrapping objects.
  class IndirectionBase
  {
    protected:
      //! Holds the type of the indirected object.
      struct t_Object
      { 
        //! Type of the indirected object.
        typedef void other; 
        //! Type of the indirected object.
        typedef const void const_other; 
      };

    public:
      //! Destructor.
      virtual ~IndirectionBase() = 0;
  
      //! Returns a reference to self.
      virtual IndirectionBase* self_() = 0; 
      //! Returns a constant reference to self.
      virtual const IndirectionBase* const self_() const = 0;
      //! Returns an object holding the type of the indirected object.
      t_Object type() const { return t_Object(); }
  };

  template < class T_OBJECT >
  class Indirection : public IndirectionBase, public T_OBJECT
  {
    protected:
      //! Holds the type of the indirected object.
      struct t_ObjectD : public IndirectionBase::t_Object
      { 
        //! Type of the indirected object.
        typedef T_OBJECT other; 
        //! Type of the indirected object.
        typedef const T_OBJECT const_other; 
      };

    public:
      //! Constructor.
      Indirection() {}
      //! Copy Constructor.
      Indirection   ( const Indirection& _c ) 
                  : T_OBJECT( _c ) {}
      //! Destructor.
      virtual ~Indirection() {}
    
      //! Returns a reference to self.
      T_OBJECT* self() { return static_cast<T_OBJECT*>(self_()); } 
      //! Returns a constant reference to self.
      const T_OBJECT* self() const { return static_cast<const T_OBJECT*>(self_()); } 
      //! Returns an object holding the type of the indirected object.
      t_ObjectD type() const { return t_ObjectD(); }
      //! Returns a reference to self.
      virtual Indirection<T_OBJECT>* self_() { return this; } 
      //! Returns a constant reference to self.
      virtual const Indirection<T_OBJECT>* self_() const 
        { return this; } 
  };
}
#endif
