//
//  Version: $Id$
//
#ifndef _TRAITS_H_
#define _TRAITS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <list>
#include <vector>
#include <string>
#include <sstream>

#include "opt/types.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

//! \brief Traits and functions capable of discerning between a few standard
//! %types and containers
//! \details The object of the Traits namespace is to offfer the ability to
//! distinguish different types when using templates. This is mostly used in
//! module \ref Genetic, where we want to write templated classes capable of
//! dealing with scalar and vectorial quantities.
namespace Traits 
{
  //! \brief Determins wether \a IS_SCALAR is indeed a scalar or vectorial property
  //! \details Declares two quantities. Dim::is_scalar and Dim::is_vector,
  //! which are respectively true and false if \a IS_SCALAR is a types::t_int,
  //! types::t_unsigned, types::t_real, char, or a bool. Otherwise,  Dim::is_scalar
  //! and Dim::is_vector are respectively false and true.
  template< class IS_SCALAR >
   struct Dim { const static bool is_scalar = false; //!< true is \a IS_SCALAR is a scalar
                const static bool is_vector = true;  //!< false if \a IS_SCALAR is a vector
   };
  //! Specialize version of Dim for types::t_real
  template<>
   struct Dim<types::t_real> { const static bool is_scalar = true;  //!< \a IS_SCALAR is a scalar
                               const static bool is_vector = false; //!< \a IS_SCALAR is not a vector
  };
  //! Specialize version of Dim for types::t_unsigned
  template<>
   struct Dim<types::t_int> { const static bool is_scalar = true;  //!< \a IS_SCALAR is a scalar
                              const static bool is_vector = false; //!< \a IS_SCALAR is not a vector
  };
  //! Specialize version of Dim for types::t_int
  template<>
   struct Dim<types::t_unsigned>
  { 
    const static bool is_scalar = true;  //!< \a IS_SCALAR is a scalar
    const static bool is_vector = false; //!< \a IS_SCALAR is not a vector
  };
  //! Specialize version of Dim for char
  template<>
   struct Dim<char> { const static bool is_scalar = true;  //!< \a IS_SCALAR is a scalar
                      const static bool is_vector = false; //!< \a IS_SCALAR is not a vector
   };
  //! Specialize version of Dim for bool
  template<>
   struct Dim<bool> { const static bool is_scalar = true;   //!< \a IS_SCALAR is a scalar
                      const static bool is_vector = false;  //!< \a IS_SCALAR is not a vector
   };

} // namespace Traits


//! \brief Templates for handling modifiers (&, *, const)
  namespace Modifier
  {
    //! Construct for detecting (absence of) const modifier
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
    //! Construct for detecting (absence of) & modifier
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
    //! Construct for detecting (absence of) * modifier
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
//
//   //! \brief Reference wrapper which forces a template algorithm to take a
//   //!        reference rather than a  value.
//   //! \details Same as reference_wrapper from boost, more or less
//   template< class T_QUANTITY >
//   class Ref
//   {
//     public:
//       typedef T_QUANTITY t_Quantity; //!< The type to reference
//
//     protected:
//       //! \cond
//       typedef typename Reference<t_Quantity> :: t_refd t_refd;
//       typedef typename Reference<t_Quantity> :: t_refd t_nonrefd;
//       //! \endcond
//
//     protected:
//       t_refd refd; //!< Reference to pass
//
//
//     public:
//       //! Constructor
//       explicit Ref( t_refd _refd ) :  refd( _refd ) {}
//
//       //! Operator returning the reference
//       t_refd operator()() { return refd; } 
//   };


    template< class T_QUANTITY > 
      typename Reference<typename Pointer<T_QUANTITY> :: t_innermost> :: t_refd
        inline innermost( T_QUANTITY &_ptr )
          { return Pointer<T_QUANTITY> :: _innermost(_ptr); }
    
    template< class T_QUANTITY > 
      typename Reference<typename Pointer<T_QUANTITY> :: t_innermost> :: t_refd
        inline const_innermost( T_QUANTITY &_ptr )
          { return Pointer<const T_QUANTITY> :: _innermost(_ptr); }
    
  }


namespace opt 
{
  //! \brief Make a vector form \a T_ARG if \a MAKEVECTOR is true.
  //! \details When setting \a MAKEVECTOR to Dim<T_ARG> :: is_vector, this
  //! function allows us to create or not a vector of T_ARG, or simply redeclare
  //! T_ARG. \relates Function::t_GradientTraits
  template< class T_ARG, bool MAKEVECTOR = true >
   struct MakeVector { typedef std::vector<T_ARG> t_Vector; //!< The the resulting type
  };
  //! Specialized version of MakeVector when \a MAKEVECTOR is set to false
  template< class T_ARG >
   struct MakeVector<T_ARG,false> { typedef T_ARG t_Vector; //!< The the resulting type
   };

  //! \brief Defines GetScalar::t_Scalar as \a T_ARG if \a T_ARG is a scalar,
  //!        and to \a T_ARG :: value_type is it is a vector.
  template< class T_ARG, bool ISVECTOR = Traits::Dim<T_ARG>::is_vector >
   struct GetScalar { typedef typename T_ARG::value_type t_Scalar; //!< the resulting type
  };
  //! Specialized version of GetScalar when ISVECTOR is set to false
  template< class T_ARG >
   struct GetScalar<T_ARG, false> { typedef T_ARG t_Scalar;  //!< The the resulting type
   };
  //! \brief Defines fuzzy math ordering operator for \a T_ARG
  //! \details if \a T_ARG is an integer type, signed or unsigned, then
  //! Fuzzy::less, Fuzzy::greater, Fuzzy::equal are exactly equivalent to <, >,
  //! and ==. If on the other hand \a T_ARG is types::real, then Fuzzy::less,
  //! Fuzzy::greater, Fuzzy::equal are defined with a fuzzy limit (types::tolerance)
  //! \param T_ARG is expected to be a scalar type
  template< class T_ARG>
    struct Fuzzy {
      //! \brief returns true if \f$ \_a  < \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or\f$ \_a  < \_b \f$
      static bool leq( const T_ARG _a, const T_ARG _b ) 
        { return _a <= _b; }
      //! \brief returns true if \f$ \_a  > \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or \f$ \_a  > \_b \f$
      static bool geq( const T_ARG _a, const T_ARG _b ) 
        { return _a >= _b; }
      //! \brief returns true if \f$ \_a  < \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  < \_b \f$
      static bool less( const T_ARG _a, const T_ARG _b ) 
        { return _a < _b; }
      //! \brief returns true if \f$ \_a  > \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
      static bool greater( const T_ARG _a, const T_ARG _b ) 
        { return _a > _b; }
      //! \brief returns true if \f$ \_a  == \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
      static bool equal( const T_ARG _a, const T_ARG _b ) 
        { return _a == _b; }
    };
  //! Specialize version for \a T_ARG a types::real, eg when fuzzyness is indeed needed 
  template<>
    struct Fuzzy<types::t_real> {
      //! \brief returns true if \f$ \_a  < \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or\f$ \_a  < \_b \f$
      static bool leq( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) < types::tolerance or _a < _b; }
      //! \brief returns true if \f$ \_a  > \_b \f$
      //! \details if \a T_ARG is a types::real, return true if 
      //! \f$|\_a - \_b| < \f$ types::tolerance or \f$ \_a  > \_b \f$
      static bool geq( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) < types::tolerance or _a > _b; }
      //! \brief returns return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  < \_b \f$
      static bool less( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) > types::tolerance and _a < _b; }
      //! \brief  return true if 
      //! \f$|\_a - \_b| > \f$ types::tolerance, and \f$ \_a  > \_b \f$
      static bool greater( const types::t_real _a, const types::t_real _b ) 
        { return std::abs(_a - _b) > types::tolerance and _a > _b; }
      //! \brief return true if \f$|\_a - \_b| > \f$ types::tolerance
      static bool equal( const types::t_real _a, const types::t_real _b ) 
        { return std::abs( _a - _b ) <  types::tolerance; }
    };

  //! \brief Defines a print_out and a broadcast function depending on 
  //! whether \a IS_SCALAR is true or false
  template<bool IS_SCALAR>
  struct IsAScalar
  {
    //! Prints out the full vector into the stream
    template< class t_Quantity >
    static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
    {
      typename t_Quantity :: const_iterator i_scal = _quantity.begin();
      typename t_Quantity :: const_iterator i_end = _quantity.end();
      for(; i_scal != i_end; ++i_scal )
        _stream << *i_scal << " ";
    }
    //! Returns a string in which is printed \a _quantity 
    template< class t_Quantity >
    static std::string print( const t_Quantity &_quantity )
    {
      std::ostringstream sstr;
      typename t_Quantity :: const_iterator i_scal = _quantity.begin();
      typename t_Quantity :: const_iterator i_end = _quantity.end();
      for(; i_scal != i_end; ++i_scal )
        sstr << *i_scal << " ";
      return sstr.str();
    }
#ifdef _MPI
    //! Broadcasts  a t_Quantity type \a _q
    template< class t_Quantity >
    static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
      { return _bc.serialize_container( _q ); }
#endif
  };
  //! Specialization of IsAScalar of \a IS_SCALAR = true
  template<>
  struct IsAScalar<true>
  {
    //! Prints out the scalar into the stream
    template< class t_Quantity >
    static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
      { _stream << _quantity; }
    //! Returns a string in which is printed \a _quantity 
    template< class t_Quantity >
    static std::string print( const t_Quantity &_quantity )
      { std::ostringstream sstr; sstr << _quantity; return sstr.str(); }
#ifdef _MPI
    //! Broadcasts  the scalar type \a _q
    template< class t_Quantity >
    static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
      { return _bc.serialize(_q); }
#endif
  };

} // namspace opt


namespace Traits
{
  //! \brief Defines a Quantity from \a T_QUANTITY, 
  //! as well as the closes scalar quantity to Quantity
  //! \details The object of the traits class is to hold for \a T_QUANTITY
  //! type, its related scalar uantity, whether it is vectorial
  //! (Quantity::is_vector) or scalar (Quantity::is_vector), as well as a number
  //! of functions which will differ if \a T_QUANTITY is truly  a scalar, or
  //! truly a vector. 
  template<class T_QUANTITY, bool ISVECTOR = Dim<T_QUANTITY> :: is_vector >
    struct Quantity  
    {
      typedef T_QUANTITY  t_Quantity;   //!< type on which to act
      //! \brief Scalar of this type
      //! \details Same as Quantity::t_Quantity if Quantity::t_Quantity is already a scalar.
      typedef typename opt::GetScalar<t_Quantity> :: t_Scalar t_ScalarQuantity; 
      //! Traits of Quantity::t_Quantity's scalar
      typedef Quantity<t_ScalarQuantity>  t_ScalarQuantityTraits;  
      //! true is Quantity::t_Quantity is a scalar 
      const static bool is_scalar = Dim<t_Quantity> :: is_scalar;
      //! true is Quantity::t_Quantity is a vector 
      const static bool is_vector = Dim<t_Quantity> :: is_vector;

      //! \brief Gets \a _n(th) component of _q 
      //! \details returns \a _q if _q is a scalar
      static t_ScalarQuantity& scalar( t_Quantity& _q, types::t_unsigned _n )
        { return _q[_n]; }
      //! \brief Gets \a _n(th) component of _q 
      //! \details returns \a _q if _q is a scalar
      static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned _n )
        { return _q[_n]; }
      //! \brief returns the size of \a _q
      //! \details returns 0 if \a _q is a scalar
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return _q.size(); }
      //! resizes \a _q, if \a is a vector
      static bool resize( const t_Quantity& _q, types::t_unsigned n ) 
        { return _q.resize(n); }
      //! Incorporates Fuzzy::less
      static bool less( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return opt::Fuzzy<t_ScalarQuantity>::less(_a,_b); }
      //! Incorporates Fuzzy::greater
      static bool greater( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return opt::Fuzzy<t_ScalarQuantity>::greater(_a,_b); }
      //! Incorporates Fuzzy::equal
      static bool equal( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return opt::Fuzzy<t_ScalarQuantity>::equal(_a,_b); }
      //! Prints out the full vector into the stream
      static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
        { opt::IsAScalar<is_scalar>::print_out(_stream, _quantity); }
      //! Prints out the full vector into a string
      static std::string print( const t_Quantity &_quantity )
        { return opt::IsAScalar<is_scalar>::print(_quantity); }
#ifdef _MPI
      //! Serializes quantity \sa mpi::BroadCast::serialize
      static bool broadcast( t_Quantity& _q, mpi::BroadCast &_bc )
        { return opt::IsAScalar<is_scalar>::broadcast(_q, _bc); }
#endif
    };

  //! \brief General traits for any function
  //! \details Defines the type of arguments and the type of the return
  template< class T_ARG, class T_RET = T_ARG >
  struct Function
  {
    typedef typename Modifier::Reference<T_ARG> :: t_nonrefd  t_Argument; //!< The Argurment type
    typedef typename Modifier::Reference<T_ARG> :: t_nonrefd  t_Return; //!< The Return type
    typedef Quantity<T_ARG> t_ArgTraits; //!< The Argument traits
    typedef Quantity<T_RET> t_RetTraits; //!< The return traits
    //! Defines a complete type for the gradient of this function
    typedef Quantity< typename opt::MakeVector< t_Return,
                          Dim<T_ARG>::is_vector > :: t_Vector > t_GradientTraits;
  };

}
#endif
