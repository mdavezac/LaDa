//
//  Version: $Id$
//
#ifndef _SEPARABLE_H_
#define _SEPARABLE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<boost/lambda/lambda.hpp>
#include<boost//boost/type_traits/function_traits.hpp>

//! Contains all things separable functions.
namespace Separable
{
  // Forward declarations
  //! \cond
  namespace details
  { 

    template< class T_BASE >
      typename T_RETIT::value_type gradient( const T_BASE &_base,
                                             typename T_BASE :: t_Arg _arg );
    template< class T_BASE, class T_ARGSIT, class T_RETIT >
      typename T_RETIT::value_type gradient( const T_BASE &_base,
                                             T_ARGSIT _arg,
                                             T_RETIT _ret  );
    //! \brief Do not use this class. Is here for type identification within
    //!        templates.
    class Base {};
    template< class T_FUNCTION,
              bool arity = boost::function_traits<T_FUNCTION> :: arity >
       class FreeFunction;

    template< class T_TYPE, bool d = false >
      T_TYPE& deref( T_TYPE );
    template< class T_TYPE, bool d = false >
      T_TYPE& deref( T_TYPE );
  }
  //! \endcond

  //! \brief Defines a separable function of many-variables.
  //! \details Programmatically the law for which this function is "separable"
  //!          is defined via \a T_GROUPOP. \a T_BASIS defines a family of 1d functions.
  //!          This function is both one-dimensional when invoked with a
  //!          scalar, and n-dimensional when invoked with iterators.
  //! \tparam T_BASIS is a container of 1d functions. These functions should
  //!         zero order evaluation via a functor call, and grdient evaluation
  //!         via a t_Return gradient( t_Arg ) member function.
  //! \tparam T_GROUPOP is template class defining how to link return values from
  //!         different basis functions together. It will be, generally,
  //!         either std::sum, or std::multiplies.
  //! \tparam T_SCALAROP is template class defining how to link return values from
  //!         a basis function with a scalar coefficient.
  template< class T_BASIS, 
  template< class T_BASIS, 
            template<class> class T_GROUPOP  = std::sum, 
            template<class> class T_SCALAROP = std::multiplies >
  class Base : public details :: Base
  {
    public:
      //! Type of the basis
      typedef T_BASIS t_Basis;
      //! Type of the operator linking a basis function and its coefficient.
      typedef T_SCALAROP< typename t_Basis::t_Return > t_ScalarOp;
      //! Type of the operator linking to basis functions.
      typedef T_GROUOP< typename t_Basis::t_Return > t_GroupOp;

    public:
      //! Type of the arguments to the one-dimensional functions.
      typedef typename t_Basis::value_type :: t_Arg t_Arg;
      //! Type of the return of the one-dimensional functions.
      typedef typename t_Basis::value_type :: t_Return t_Return;

    public:
      //! Constructor
      Base() : basis() { coefs.resize( basis.size() ); }
      //! Destructor
      ~Base() {}

      //! Returns the function evaluated at \a _args.
      template< template<class> class T_CONTAINER >
        t_Return operator()( const T_CONTAINER<t_Args> &_args ) const
        { return operator()( _args.begin(), _args.end() ); }
      //! \brief Returns the function evaluated by variables in range 
      //!        \a _first to \a _last.
      //! \details In this case, the function is of as many variables as there
      //!          are functions in the basis.
      template< class T_ITERATOR >
        t_Return operator()( T_ITERATOR _first, T_ITERATOR _last ) const;
      //! \brief Returns the gradient of the one dimensional function.
      t_Return gradient( t_Arg _arg ) const 
        { return details::gradient( *this, _arg ); }
      //! Computes the gradient and stores it in \a _ret.
      template< class T_ARGIT, class T_RETIT >
        void gradient( T_ARGIT _first, T_RETIT _ret ) const
        { return details::gradient( *this, _first, _ret ); }
      //! \brief Return the function evaluated at \a _arg
      t_Return operator()( t_Arg _arg ) const;
      //! Returns a reference to coeficient _i
      t_Return& operator[]( types::t_unsigned _i ) { return coefs[_i]; }
      //! Returns a reference to coeficient _i
      const t_Return& operator[]( types::t_unsigned _i ) const 
        { return coefs[_i]; }
      //! Sets all coefs in range \a _first to \a _last.
      template< class T_ITERATOR >
        void set( T_ITERATOR _first, T_ITERATOR _last )
      //! Returns a reference to the basis
      t_Basis& Basis() { return basis; }
      //! Returns a constant reference to the basis.
      const t_Basis& Basis() const { return basis; }
      //! Serializes a structure.
      template<class ARCHIVE> void serialize( ARCHIVE & _ar,
                                              const unsigned int _version);
     
    protected:
      //! A family of functions. 
      Basis basis;
      //! \brief Type of the contained for the coefficients of a single
      //!        separable function.
      typedef std::vector< t_Return > t_Coefs;
      //! A container of coefficients.
      std::vector< t_Return > coefs;
      //! Links basis functions.
      t_GroupOp groupop;
      //! Links scalars to basis functions.
      t_ScalarOp scalarop;

    protected:
      //! For alternating least-square purposes.
      template< class T_ITIN, class T_ITOUT>
        t_Return expand( T_ITIN _first, T_ITIN _last, T_ITOUT _out ) const;
  };

  // Forward declaration.
  template< class T_ALLSQ > class AllsqInterface;

  //! Factor of a separable function
  //! \tparam T_BASIS is a container of 1d functions. These functions should
  //!         zero order evaluation via a functor call, and grdient evaluation
  //!         via a t_Return gradient( t_Arg ) member function.
  template< class T_BASIS > class Factor : public Base< T_BASIS >{};
  //! One single separable function.
  //! \tparam T_BASIS is a container of 1d functions. These functions should
  //!         zero order evaluation via a functor call, and grdient evaluation
  //!         via a t_Return gradient( t_Arg ) member function.
  template< class T_BASIS > class Summand :
      public Base< std::vector< Factor<T_BASIS, std::multiplies, std::sum> > >{};
  //! \brief A sum of separable functions.
  //! \tparam T_BASIS is a container of 1d functions. These functions should
  //!         zero order evaluation via a functor call, and grdient evaluation
  //!         via a t_Return gradient( t_Arg ) member function.
  template< class T_BASIS >
    class Function : public Base< std::vector< Summand<T_BASIS> > >
    {
        template<class T_ALLSQ>  friend class AllsqInterface;
        //! Type of the base class.
        typedef Base< std::vector< Summand<T_BASIS> > > t_Base;
      public:
        //! Type of the basis
        typedef T_BASIS t_Basis;
        //! Type of the arguments to the one-dimensional functions.
        typedef typename t_Basis::value_type :: t_Arg t_Arg;
        //! Type of the return of the one-dimensional functions.
        typedef typename t_Basis::value_type :: t_Return t_Return;
 
      public:
        //! Constructor
        Function() : t_Base() {}
        //! Destructor
        ~Function() {}
 
        //! Returns the function evaluated at \a _args.
        template< template<class> class T_CONTAINER >
          t_Return operator()( const T_CONTAINER<t_Args> &_args ) const
          { return operator()( _args.begin(), _args.end() ); }
        //! \brief Returns the function evaluated by variables in range 
        //!        \a _first to \a _last.
        //! \details In this case, the function is of as many variables as there
        //!          are functions in the basis.
        template< class T_ITERATOR >
          t_Return operator()( T_ITERATOR _first, T_ITERATOR _last ) const;
        //! Computes the gradient and stores in \a _ret.
        template< class T_ARGIT, class T_RETIT >
          void gradient( T_ARGIT _first, T_RETIT _ret ) const;
 
      protected:
        using t_Base::groupop;
        using t_Base::scalarop;
    };

  //! \brief Creates the collapsed matrices and vectors for the
  //!        alternating linear least-square fit method.
  //! \param[in] _func a function for which to map out the Altenating
  //!             linear least-square fit matrices
  //! \param[out] The Alternating linear least-square fit matrix.
  //!             It is organized as follows:
  //!             Each dimension dependent sub-matrix is one value in the vector _A.
  //!             each sub-matrix is a vector organized as follows
  //!               obs 0     obs 1
  //!             yyyyyyyy  yyyyyyyyyy
  //!             Each yyyyyyy runs as 
  //!               r=0    r=1  rank of separable function
  //!              xxxxx  xxxxx 
  //!             Each xxxxx as 
  //!               m=0     m=1   basis for factor n of separable function.
  //! \tparam T_BASIS is the 1d basis for the separable functions.
  template< class T_ALLSQ, class T_BASIS >
    void assign_Allsq_matrix( const Function<T_BASIS> &_func,
                              const T_ALLSQ::t_Vectors &_input,
                             T_ALLSQ::t_Matrices &_A )
  //! \brief assigns coefficients computed from an Alternating linear
  //         least-square fit to \a _func.
  template<class T_ALLSQ, class T_BASIS >
    void assign_from_allsq( Function<T_BASIS> &_func, const t_Vectors &_coefs );

}

#include "separable.impl.h"

#endif
