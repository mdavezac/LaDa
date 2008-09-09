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
        typedef void type; 
      }

    public:
      //! Constructor.
      virtual IndirectionBase() = 0;
      //! Destructor.
      virtual ~IndirectionBase() = 0;
  
      //! Returns a reference to self.
      virtual IndirectionBase* self() = 0; 
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
      struct t_ObjectD
      { 
        //! Type of the indirected object.
        typedef T_OBJECT type; 
      }

    public:
      //! Constructor.
      virtual Indirection() {}
      //! Copy Constructor.
      virtual Indirection   ( const Indirection& _c ) 
                          : T_OBJECT( _c ) {}
      //! Destructor.
      virtual ~Indirection() {}
    
      //! Returns a reference to self.
      Indirection* self() { return static_cast< T_OBJECT*>( this ); } 
      //! Returns a constant reference to self.
      virtual const Indirection* self() const 
        { return static_cast< const T_OBJECT* const >( this ); } 
      //! Returns an object holding the type of the indirected object.
      t_ObjectD type() const { return t_ObjectD(); }
  };
}

