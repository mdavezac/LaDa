#ifndef _OPT_MODIFIERS_H_
#define _OPT_MODIFIERS_H_

#include "LaDaConfig.h"

namespace LaDa
{
  //! Templates for handling modifiers (&, *, const)
  namespace Modifier
  {
    //! Construct for detecting (absence/presence of) const modifier
    template<class T_QUANTITY> struct Const
    { 
      const static bool is_const = false;  //!< True ie T_QUANTITY is const
      typedef T_QUANTITY t_nonconstant;    //!< Non Constant type
      typedef const T_QUANTITY t_constant; //!< Constant type
    };
    //! Construct for detecting (presence of) const modifier
    template<class T_QUANTITY> struct Const<const T_QUANTITY>
    { 
      const static bool is_const = true;  //!< True ie T_QUANTITY is const
      typedef T_QUANTITY t_nonconstant;    //!< Non Constant type
      typedef const T_QUANTITY t_constant; //!< Constant type
    };
    //! Construct for detecting (absence/presence of) & modifier
    template<class T_QUANTITY> struct Reference
    { 
      const static bool is_refd = false;  //!< True ie T_QUANTITY is a reference
      typedef T_QUANTITY t_nonrefd;    //!< Non Ref'd type
      typedef T_QUANTITY& t_refd;      //!< Ref'd type
    };
    //! Construct for detecting (presence of) & modifier
    template<class T_QUANTITY> struct Reference<T_QUANTITY&>
    { 
      const static bool is_refd = true;  //!< True ie T_QUANTITY is const
      typedef T_QUANTITY t_nonrefd;      //!< Non ref'd type
      typedef T_QUANTITY& t_refd;        //!< Ref'd type
    };
    //! Construct for detecting (absence/presence of) * modifier
    template<class T_QUANTITY> class Pointer
    { 
      //! \cond
      typedef typename Reference<T_QUANTITY> :: t_refd t_by_address;
      template< class TQUANTITY > 
        typename Reference<typename Pointer<TQUANTITY> :: t_innermost> :: t_refd
          inline innermost( TQUANTITY &_ptr );
      //! \endcond

      public:
        typedef T_QUANTITY t_Quantity; //!< Original type

      public:
        const static bool is_pointer = false;  //!< True ie T_QUANTITY is a reference
        typedef T_QUANTITY t_nonpointer;       //!< Non-pointer type
        typedef T_QUANTITY* t_pointer;         //!< Pointer type
        typedef T_QUANTITY t_innermost;  //!< Innermost pointerd type

      public:
        //! Returns reference to innermost value of pointer
        static typename Reference<t_nonpointer>::t_refd 
          _innermost( t_by_address _ptr ) { return _ptr; }
    };
    //! Construct for detecting (presence of) * modifier
    template<class T_QUANTITY> class Pointer<T_QUANTITY*>
    { 
      //! \cond
      typedef T_QUANTITY* const t_by_address;
      template< class TQUANTITY > 
        typename Reference<typename Pointer<TQUANTITY> :: t_innermost> :: t_refd
          inline innermost( TQUANTITY &_ptr );
      //! \endcond

      public:
        typedef T_QUANTITY t_Quantity; //!< Original type

      public:
        const static bool is_pointer = true;  //!< True ie T_QUANTITY is a reference
        typedef T_QUANTITY t_nonpointer;       //!< Non-pointer type
        typedef T_QUANTITY* t_pointer;         //!< Pointer type
         //!  Innermost pointed type
        typedef typename Pointer<t_nonpointer>::t_innermost t_innermost; 

      public:
        //! Returns reference to innermost value of pointer
        static typename Reference< t_innermost >::t_refd 
          _innermost( t_by_address _ptr )
          { return Pointer<t_nonpointer>::_innermost(*_ptr); }
    };
    //! Construct for detecting (presence of) * const modifier
    template<class T_QUANTITY> class Pointer<const T_QUANTITY>
    { 
      //! \cond
      typedef const T_QUANTITY& t_by_address;
      template< class TQUANTITY > 
        typename Reference<typename Pointer<TQUANTITY> :: t_innermost> :: t_refd
          inline innermost( TQUANTITY &_ptr );
      //! \endcond

      public:
        typedef T_QUANTITY t_Quantity; //!< Original type

      public:
        const static bool is_pointer = true;  //!< True ie T_QUANTITY is a reference
        typedef T_QUANTITY const t_nonpointer;       //!< Non-pointer type
        typedef T_QUANTITY* const t_pointer;         //!< Pointer type
         //!  Innermost pointed type
        typedef T_QUANTITY const t_innermost; 

      public:
        //! Returns reference to innermost value of pointer
        static typename Reference< t_nonpointer >::t_refd 
          _innermost( t_by_address _ptr ) { return _ptr; }
    };
    //! Construct for detecting (presence of) * const modifier
    template<class T_QUANTITY> class Pointer<T_QUANTITY* const>
    { 
      //! \cond
      typedef T_QUANTITY* t_by_address;
      template< class TQUANTITY > 
        typename Reference<typename Pointer<TQUANTITY> :: t_innermost> :: t_refd
          inline innermost( TQUANTITY &_ptr );
      //! \endcond

      public:
        typedef T_QUANTITY t_Quantity; //!< Original type

      public:
        const static bool is_pointer = true;  //!< True ie T_QUANTITY is a reference
        typedef T_QUANTITY const t_nonpointer;       //!< Non-pointer type
        typedef T_QUANTITY* const t_pointer;         //!< Pointer type
         //!  Innermost pointed type
        typedef typename Pointer<t_nonpointer>::t_innermost t_innermost; 

      public:
        //! Returns reference to innermost value of pointer
        static typename Reference< t_innermost >::t_refd 
          _innermost( t_by_address _ptr )
          { return Pointer<t_nonpointer>::_innermost(*_ptr); }
    };

    //! \brief Helper %function returning the value to which a pointer points.
    //! \details Say you call this %function with an object in argument, you get
    //!          a reference to this object in return. Say you call it with a
    //!          pointer to the same object. You still get the same reference in
    //!          return. Now if you call it with a pointer to this pointer to
    //!          the same object, then  in that case, you still get the same
    //!          reference. In other words. whatever you put in you get the most
    //!          dereferenced value. Doesn't dereference iterators though...
    template< class T_QUANTITY > 
      typename Reference<typename Pointer<T_QUANTITY> :: t_innermost> :: t_refd
        inline innermost( T_QUANTITY &_ptr )
          { return Pointer<T_QUANTITY> :: _innermost(_ptr); }
    
    //! \brief const version of Modifier::innnermost()
    template< class T_QUANTITY > 
      typename Reference<const typename Pointer<T_QUANTITY> :: t_innermost> :: t_refd
        inline const_innermost(const T_QUANTITY &_ptr )
          { return Pointer<const T_QUANTITY> :: _innermost(_ptr); }

    //! If then else. True flavor.
    template< bool CONDITION, class THEN, class ELSE >
      struct if_then_else 
      {
        //! result.
        typedef THEN type;
      };
    template< class THEN, class ELSE >
      struct if_then_else< false, THEN, ELSE >
      {
        //! result.
        typedef ELSE type;
      };
    
  }
} // namespace LaDa
#endif
