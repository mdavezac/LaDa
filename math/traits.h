#ifndef _OPT_TRAITS_H_
#define _OPT_TRAITS_H_

#include "LaDaConfig.h"


#include <list>
#include <vector>
#include <string>
#include <sstream>

#include <mpi/mpi_object.h>
#include <opt/types.h>

#include "fuzzy.h"
#include "modifiers.h"


namespace LaDa
{
  namespace math
  {
    //! \brief %Traits and functions capable of discerning between a few standard
    //! %types and containers
    //! \details The object of the Traits namespace is to offfer the ability to
    //! distinguish different %types when using templates. This is mostly used in
    //! module \ref Genetic, where we want to write templated classes capable of
    //! dealing with scalar and vectorial quantities.
    namespace traits 
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
      //! Specialize version of Dim for types::t_complex
      template<>
       struct Dim< types::t_complex >
       { 
         const static bool is_scalar = true;   //!< \a IS_SCALAR is a scalar
         const static bool is_vector = false;  //!< \a IS_SCALAR is not a vector
       };
      //! Specialize version of Dim for std::string
      template<>
       struct Dim< std::string >
       { 
         const static bool is_scalar = true;   //!< \a IS_SCALAR is a scalar
         const static bool is_vector = false;  //!< \a IS_SCALAR is not a vector
       };
      //! Specialize version of Dim for constant %types
      template<class T_QUANTITY>
       struct Dim<const T_QUANTITY>
       {
         const static bool is_scalar = Dim<T_QUANTITY>::is_scalar;  //!< \a IS_SCALAR is a scalar
         const static bool is_vector = Dim<T_QUANTITY>::is_vector; //!< \a IS_SCALAR is not a vector
       };
 
    } // namespace Traits
 
 
    //! \brief Make a vector form \a T_ARG if \a MAKEVECTOR is true.
    //! \details When setting \a MAKEVECTOR to Dim<T_ARG> :: is_vector, this
    //! %function allows us to create or not a vector of T_ARG, or simply redeclare
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
    template< class T_ARG, bool ISVECTOR = traits::Dim<T_ARG>::is_vector >
     struct GetScalar { typedef typename T_ARG::value_type t_Scalar; //!< the resulting type
    };
    //! Specialized version of GetScalar when ISVECTOR is set to false
    template< class T_ARG >
     struct GetScalar<T_ARG, false> { typedef T_ARG t_Scalar;  //!< The the resulting type
     };
 
 
    //! \brief Defines a print_out and a serialize %function depending on 
    //! whether \a IS_SCALAR is true or false
    template<bool IS_SCALAR>
    struct IsAScalar
    {
      //! Prints out the full vector into the stream
      template< class t_Quantity >
      static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
      {
        typedef traits::Dim<typename t_Quantity :: value_type> t_Dim;
        typedef IsAScalar< t_Dim :: is_scalar > t_IsAScalar;
        typename t_Quantity :: const_iterator i_scal = _quantity.begin();
        typename t_Quantity :: const_iterator i_end = _quantity.end();
        for(; i_scal != i_end; ++i_scal )
          _stream <<  t_IsAScalar :: print(*i_scal) << " ";
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
      //! \brief Gets \a _n(th) component of _q 
      //! \details returns \a _q if _q is a scalar
      template< class t_Quantity >
      static t_Quantity& scalar( t_Quantity& _q, types::t_unsigned _n )
        { return _q[_n]; }
      //! \brief Gets \a _n(th) component of _q 
      //! \details returns \a _q if _q is a scalar
      template< class t_Quantity >
      static const t_Quantity& scalar( const t_Quantity& _q, types::t_unsigned _n )
        { return _q[_n]; }
      //! \brief returns the size of \a _q
      //! \details returns 0 if \a _q is a scalar
      template< class t_Quantity >
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return _q.size(); }
      //! resizes \a _q, if \a is a vector
      template< class t_Quantity >
      static bool resize( const t_Quantity& _q, types::t_unsigned n ) 
        { return _q.resize(n); }
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
      //! \brief Gets \a _n(th) component of _q 
      //! \details returns \a _q if _q is a scalar
      template< class t_Quantity >
      static t_Quantity& scalar( t_Quantity& _q, types::t_unsigned _n )
        { return _q; }
      //! \brief Gets \a _n(th) component of _q 
      //! \details returns \a _q if _q is a scalar
      template< class t_Quantity >
      static const t_Quantity& scalar( const t_Quantity& _q, types::t_unsigned _n )
        { return _q; }
      //! \brief returns the size of \a _q
      //! \details returns 1 if \a _q is a scalar
      template< class t_Quantity >
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return 1; }
      //! resizes \a _q, if \a is a vector
      template< class t_Quantity >
      static bool resize( const t_Quantity& _q, types::t_unsigned n ) 
        { return true; }
    };
 
    namespace traits
    {
      //! \brief Defines %types pertaining to \a T_QUANTITY, 
      //! eg itself, its components, ...
      //! \details The object of the traits class is to hold for \a T_QUANTITY
      //! type, its related scalar uantity, whether it is vectorial
      //! (Quantity::is_vector) or scalar (Quantity::is_vector), as well as a number
      //! of functions which will differ if \a T_QUANTITY is truly  a scalar, or
      //! truly a vector. 
      template<class T_QUANTITY, bool ISVECTOR = Dim<T_QUANTITY> :: is_vector >
        struct Quantity  
        {
          typedef T_QUANTITY  t_Quantity;   //!< type on which to act
          //! constant quantity type
          typedef typename Modifier::Const<t_Quantity> :: t_constant t_const_Quantity;
          //! \brief Scalar of this type
          //! \details Same as Quantity::t_Quantity if Quantity::t_Quantity is already a scalar.
          typedef typename GetScalar<t_Quantity> :: t_Scalar t_ScalarQuantity; 
          //! traits of Quantity::t_Quantity's scalar
          typedef Quantity<t_ScalarQuantity>  t_ScalarQuantityTraits;  
          //! traits of constant Quantity::t_Quantity
          typedef Quantity< t_const_Quantity >  t_const_QuantityTraits;  
          //! true is Quantity::t_Quantity is a scalar 
          const static bool is_scalar = Dim<t_Quantity> :: is_scalar;
          //! true is Quantity::t_Quantity is a vector 
          const static bool is_vector = Dim<t_Quantity> :: is_vector;
 
          //! Incorporates math::le
          static bool le( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
            { return math::le(_a,_b); }
          //! Incorporates math::leq
          static bool leq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
            { return math::leq(_a,_b); }
          //! Incorporates math::gt
          static bool gt( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
            { return math::gt(_a,_b); }
          //! Incorporates math::geq
          static bool geq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
            { return math::geq(_a,_b); }
          //! Incorporates math::eq
          static bool eq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
            { return math::eq(_a,_b); }
          //! Incorporates math::neq
          static bool neq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
            { return not math::eq(_a,_b); }
          //! Prints out the full vector into the stream
          static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
            { IsAScalar<is_scalar>::print_out(_stream, _quantity); }
          //! Prints out the full vector into a string
          static std::string print( const t_Quantity &_quantity )
            { return IsAScalar<is_scalar>::print(_quantity); }
          //! \brief Gets \a _n(th) component of _q 
          //! \details returns \a _q if _q is a scalar
          static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned _n )
            { return IsAScalar<is_scalar>::scalar(_q, _n); }
          //! \brief returns the size of \a _q
          //! \details returns 0 if \a _q is a scalar
          static types::t_unsigned size( const t_Quantity& _q ) 
            { return IsAScalar<is_scalar>::size(_q); }
          //! resizes \a _q, if \a is a vector
          static bool resize( const t_Quantity& _q, types::t_unsigned n ) 
            { return IsAScalar<is_scalar>::resize(_q); }
        };
 
 
      //! \brief General traits for any %function
      //! \details Defines the type of arguments and the type of the return
      template< class T_ARG, class T_RET = T_ARG >
      struct Function
      {
        //! The Argurment type
        typedef typename Modifier::Reference<T_ARG> :: t_nonrefd  t_Argument; 
        //! The Return type
        typedef typename Modifier::Reference<T_ARG> :: t_nonrefd  t_Return; 
        typedef Quantity<T_ARG> t_ArgTraits; //!< The Argument traits
        typedef Quantity<T_RET> t_RetTraits; //!< The return traits
        //! Defines a complete type for the gradient of this %function
        typedef Quantity< typename MakeVector< t_Return,
                              Dim<T_ARG>::is_vector > :: t_Vector > t_GradientTraits;
      };
 
 
    } // namespace traits
  } // namespace math
} // namespace Lada

#endif
