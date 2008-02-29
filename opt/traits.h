//
//  Version: $Id$
//
#ifndef _OPT_TRAITS_H_
#define _OPT_TRAITS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <list>
#include <vector>
#include <string>
#include <sstream>

#include "types.h"
#include "fuzzy.h"
#include "modifiers.h"

#include <mpi/mpi_object.h>

//! \brief %Traits and functions capable of discerning between a few standard
//! %types and containers
//! \details The object of the Traits namespace is to offfer the ability to
//! distinguish different %types when using templates. This is mostly used in
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


namespace opt 
{
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
  template< class T_ARG, bool ISVECTOR = Traits::Dim<T_ARG>::is_vector >
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
      typedef Traits::Dim<typename t_Quantity :: value_type> t_Dim;
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
#ifdef _MPI
    /** \ingroup MPI
     * \brief Serializes  a t_Quantity type \a _q 
     **/
    template< class t_Quantity >
    static bool serialize( t_Quantity& _q, mpi::BroadCast &_bc )
      { return mpi::BroadCast::Container<t_Quantity>(_bc)( _q ); }
    /** \ingroup MPI
     * \brief Serializes  a constant t_Quantity type \a _q
     **/
    template< class t_Quantity >
    static bool serialize( const t_Quantity& _q, mpi::BroadCast &_bc )
      { return mpi::BroadCast::Container<t_Quantity>(_bc)( _q ); }
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
#ifdef _MPI
    //! Broadcasts  the scalar type \a _q
    template< class t_Quantity >
    static bool serialize( t_Quantity& _q, mpi::BroadCast &_bc )
      { return _bc.serialize(_q); }
    //! Broadcasts  the constant scalar type \a _q
    template< class t_Quantity >
    static bool serialize( const t_Quantity& _q, mpi::BroadCast &_bc )
      { return _bc.serialize(_q); }
#endif
  };

} // namspace opt


namespace Traits
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
      typedef typename opt::GetScalar<t_Quantity> :: t_Scalar t_ScalarQuantity; 
      //! Traits of Quantity::t_Quantity's scalar
      typedef Quantity<t_ScalarQuantity>  t_ScalarQuantityTraits;  
      //! Traits of constant Quantity::t_Quantity
      typedef Quantity< t_const_Quantity >  t_const_QuantityTraits;  
      //! true is Quantity::t_Quantity is a scalar 
      const static bool is_scalar = Dim<t_Quantity> :: is_scalar;
      //! true is Quantity::t_Quantity is a vector 
      const static bool is_vector = Dim<t_Quantity> :: is_vector;

      //! Incorporates Fuzzy::le
      static bool le( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Fuzzy::le(_a,_b); }
      //! Incorporates Fuzzy::leq
      static bool leq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Fuzzy::leq(_a,_b); }
      //! Incorporates Fuzzy::gt
      static bool gt( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Fuzzy::gt(_a,_b); }
      //! Incorporates Fuzzy::geq
      static bool geq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Fuzzy::geq(_a,_b); }
      //! Incorporates Fuzzy::eq
      static bool eq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return Fuzzy::eq(_a,_b); }
      //! Incorporates Fuzzy::neq
      static bool neq( const t_ScalarQuantity _a, const t_ScalarQuantity _b ) 
        { return not Fuzzy::eq(_a,_b); }
      //! Prints out the full vector into the stream
      static void print_out( std::ostream& _stream, const t_Quantity &_quantity )
        { opt::IsAScalar<is_scalar>::print_out(_stream, _quantity); }
      //! Prints out the full vector into a string
      static std::string print( const t_Quantity &_quantity )
        { return opt::IsAScalar<is_scalar>::print(_quantity); }
      //! \brief Gets \a _n(th) component of _q 
      //! \details returns \a _q if _q is a scalar
      static const t_ScalarQuantity& scalar( const t_Quantity& _q, types::t_unsigned _n )
        { return opt::IsAScalar<is_scalar>::scalar(_q, _n); }
      //! \brief returns the size of \a _q
      //! \details returns 0 if \a _q is a scalar
      static types::t_unsigned size( const t_Quantity& _q ) 
        { return opt::IsAScalar<is_scalar>::size(_q); }
      //! resizes \a _q, if \a is a vector
      static bool resize( const t_Quantity& _q, types::t_unsigned n ) 
        { return opt::IsAScalar<is_scalar>::resize(_q); }
#ifdef _MPI
      //! Serializes quantity \sa mpi::BroadCast::serialize
      static bool serialize( t_Quantity& _q, mpi::BroadCast &_bc )
        { return opt::IsAScalar<is_scalar>::serialize(_q, _bc); }
      //! Serializes constant quantity \sa mpi::BroadCast::serialize
      static bool serialize( const t_Quantity& _q, mpi::BroadCast &_bc )
        { return opt::IsAScalar<is_scalar>::serialize(_q, _bc); }
#endif
    };


  //! \brief General traits for any %function
  //! \details Defines the type of arguments and the type of the return
  template< class T_ARG, class T_RET = T_ARG >
  struct Function
  {
    typedef typename Modifier::Reference<T_ARG> :: t_nonrefd  t_Argument; //!< The Argurment type
    typedef typename Modifier::Reference<T_ARG> :: t_nonrefd  t_Return; //!< The Return type
    typedef Quantity<T_ARG> t_ArgTraits; //!< The Argument traits
    typedef Quantity<T_RET> t_RetTraits; //!< The return traits
    //! Defines a complete type for the gradient of this %function
    typedef Quantity< typename opt::MakeVector< t_Return,
                          Dim<T_ARG>::is_vector > :: t_Vector > t_GradientTraits;
  };


}


namespace mpi
{
  template< class T_QUANTITY >
  void send_quantity( types::t_unsigned _bull,
                      const T_QUANTITY &_q, MPI::Intracomm *_comm)
  {
    if( _comm->Get_rank() < 2 ) return;
    typedef Traits::Quantity< T_QUANTITY > t_QuantityTraits;
    BroadCast bc( _comm );
    t_QuantityTraits :: serialize(_q, bc );
    bc.allocate_buffers();
    t_QuantityTraits :: serialize(_q, bc );
    bc.send_ptp( _bull );
  }

  template< class T_QUANTITY >
  void receive_quantity( types::t_unsigned _bull,
                         T_QUANTITY &_q, MPI::Intracomm *_comm)
  {
    if( _comm->Get_rank() < 2 ) return;
    typedef Traits::Quantity< T_QUANTITY > t_QuantityTraits;
    BroadCast bc( _comm );
    t_QuantityTraits :: serialize(_q, bc );
    bc.allocate_buffers();
    t_QuantityTraits :: serialize(_q, bc );
    bc.receive_ptp( _bull );
    t_QuantityTraits :: serialize(_q, bc );
  }
}

#endif
