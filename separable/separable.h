//
//  Version: $Id$
//
#ifndef _SEPARABLE_H_
#define _SEPARABLE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<functional>
#include<numeric>

#include<boost/lambda/lambda.hpp>
#include<boost/type_traits/function_traits.hpp>

#include<opt/types.h>
#include<opt/debug.h>

#include "collapse.h"

//! Contains all things separable functions.
namespace Separable
{
  // Forward declarations
  //! \cond
  template< class T_FUNCTION >  class Collapse;
  template< class T_B, template<class> class T_G, template<class> class T_S > class Base;
  namespace details
  { 
    template< class T_FUNCTION,
              bool arity = boost::function_traits<T_FUNCTION> :: arity >
       class FreeFunction;
  }
  //! \endcond

// //! Prints out function.
// template< class T_B>  std::ostream& operator<<( std::ostream& _stream,
//                                                 const Function<T_B>& _func );
  //! Prints out function.
  template< class T_B, template<class> class T_G, template<class> class T_S>
  std::ostream& operator<<( std::ostream& _stream,
                            const Base<T_B, T_G, T_S>& _func );
 
  //! \brief Defines a separable function of many-variables.
  //! \details Programmatically the law for which this function is "separable"
  //!          is defined via \a T_GROUPOP. \a T_BASIS defines a family of 1d functions.
  //!          This function is both one-dimensional when invoked with a
  //!          scalar, and n-dimensional when invoked with iterators.
  //! \tparam T_BASIS is a container of 1d functions. These functions should
  //!         zero order evaluation via a functor call.
  //! \tparam T_GROUPOP is template class defining how to link return values from
  //!         different basis functions together. It will be, generally,
  //!         either std::plus, or std::multiplies.
  //! \tparam T_SCALAROP is template class defining how to link return values from
  //!         a basis function with a scalar coefficient.
  template< class T_BASIS, 
            template<class> class T_GROUPOP  = std::plus, 
            template<class> class T_SCALAROP = std::multiplies >
  class Base 
  {
      template< class T_FUNCTION > friend class Collapse;
      friend std::ostream& operator<< <T_BASIS, T_GROUPOP, T_SCALAROP>
                                      ( std::ostream& _stream,
                                        const Base< T_BASIS, 
                                                    T_GROUPOP,
                                                    T_SCALAROP >& _func );
    public:
      //! Type of the basis
      typedef T_BASIS t_Basis;
      //! Type of the arguments to the one-dimensional functions.
      typedef typename t_Basis::value_type :: t_Arg t_Arg;
      //! Type of the return of the one-dimensional functions.
      typedef typename t_Basis::value_type :: t_Return t_Return;
      //! Type of the operator linking a basis function and its coefficient.
      typedef T_SCALAROP< t_Return > t_ScalarOp;
      //! Type of the operator linking to basis functions.
      typedef T_GROUPOP< t_Return > t_GroupOp;
      //! \brief Type of the contained for the coefficients of a single
      //!        separable function.
      typedef std::vector< t_Return > t_Coefs;

      //! Whether this function has gradients
      const static bool has_gradient;
      //! A family of functions. 
      t_Basis basis;
      //! A container of coefficients.
      t_Coefs coefs;

      //! Constructor
      Base() : basis(), name("Function") { coefs.resize( basis.size() ); }
      //! Destructor
      ~Base() {}

      //! Returns the function evaluated at \a _args.
      template< template<class> class T_CONTAINER >
        t_Return operator()( const T_CONTAINER<t_Arg> &_args ) const
        { return operator()( _args.begin(), _args.end() ); }
      //! \brief Returns the function evaluated by variables in range 
      //!        \a _first to \a _last.
      //! \details In this case, the function is of as many variables as there
      //!          are functions in the basis.
      template< class T_ITERATOR >
        t_Return operator()( T_ITERATOR _first, T_ITERATOR _last ) const;
      //! \brief Return the function evaluated at \a _arg
      t_Return operator()( t_Arg _arg ) const;
      //! Returns a reference to coeficient _i
      t_Return& operator[]( types::t_unsigned _i ) { return coefs[_i]; }
      //! Returns a reference to coeficient _i
      const t_Return& operator[]( types::t_unsigned _i ) const 
        { return coefs[_i]; }
      //! Sets all coefs in range \a _first to \a _last.
      template< class T_ITERATOR >
        void set( T_ITERATOR _first, T_ITERATOR _last );
      //! Serializes a structure.
      template<class ARCHIVE>
        void serialize( ARCHIVE & _ar, const unsigned int _version);
      //! Sets name of function.
      void set_name( const std::string &_str ) { name = _str; }

    protected:
      //! Links basis functions.
      t_GroupOp groupop;
      //! Links scalars to basis functions.
      t_ScalarOp scalarop;

    protected:
      //! \brief Name of this stage of the separable function.
      //! \details Should be set by the sum of separable function.
      std::string name;
  };
  
  template< class T_BASIS, 
            template<class> class T_GROUPOP,
            template<class> class T_SCALAROP>
   const bool Base<T_BASIS, T_GROUPOP, T_SCALAROP> :: has_gradient = T_BASIS::has_gradient;

  // Forward declaration.
  template< class T_ALLSQ > class AllsqInterface;

  //! Factor of a separable function
  //! \tparam T_BASIS is a container of 1d functions. These functions should
  //!         zero order evaluation via a functor call, and grdient evaluation
  //!         via a t_Return gradient( t_Arg ) member function.
  template< class T_BASIS > class Factor : public Base< T_BASIS >{};
  //! One single separable function.
  //! \tparam T_BASIS is a container of 1d functions. These functions should
  //!         zero order evaluation via a functor call.
  template< class T_BASIS > class Summand :
      public Base< std::vector< Factor<T_BASIS> >, std::multiplies, std::plus >{};
  /** \brief A sum of separable functions.
   * \details The separable function \f$F(\mathbf{x})\f$ acting upon vector
   *          \f$\mathbf{x}\f$ and returning a scalar can be defined as 
   *    \f[
   *        F(\mathbf{x}) = \sum_r \prod_d \sum_i \lambda_{d,i}^{(r)}
   *                        g_{i,n}^{(r)}(x_i),
   *    \f]
   *          with the sum over \e r running over the ranks of the separable
   *          functions, the product over \e i are the separable functions
   *          proper, and the sum_n is an expansion of the factors ver some
   *          family of 1d-functions.
   * \tparam T_BASIS is a container of 1d functions. These functions should
   *         return a zero order evaluation via a functor call. **/
  template< class T_BASIS >
    class Function : public Base< T_BASIS >
    {
//       friend std::ostream& operator<< <T_BASIS>( std::ostream& _stream,
//                                                  const Function<T_BASIS>& _func );
        friend class Collapse< Function<T_BASIS> >;
        template<class T_ALLSQ>  friend class AllsqInterface;
        //! Type of the base class.
        typedef Base< T_BASIS > t_Base;
      public:
        //! Type of the basis
        typedef T_BASIS t_Basis;
        //! Type of the arguments to the one-dimensional functions.
        typedef typename t_Basis :: value_type :: t_Arg t_Arg;
        //! Type of the return of the one-dimensional functions.
        typedef typename t_Basis :: value_type :: t_Return t_Return;
 
      public:
        //! Constructor
        Function() : t_Base() {}
        //! Destructor
        ~Function() {}
 
        //! Returns the function evaluated at \a _args.
        template< template<class> class T_CONTAINER >
          t_Return operator()( const T_CONTAINER<t_Arg> &_args ) const
          { return operator()( _args.begin(), _args.end() ); }
        //! \brief Returns the function evaluated by variables in range 
        //!        \a _first to \a _last.
        //! \details In this case, the function is of as many variables as there
        //!          are functions in the basis.
        template< class T_ITERATOR >
          t_Return operator()( T_ITERATOR _first, T_ITERATOR _last ) const;
 
      protected:
        //! \cond 
        using t_Base :: groupop;
        using t_Base :: scalarop;
        using t_Base :: coefs;
        using t_Base :: basis;
        using t_Base :: has_gradient;
        typedef typename t_Base :: t_Coefs t_Coefs;
        //! \endconf
    };

} // end of Separable namespace

#include "separable.impl.h"

#endif
