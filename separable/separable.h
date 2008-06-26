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

//! Contains all things separable functions.
namespace Separable
{
  // Forward declarations
  //! \cond
  template< class T_FUNCTION >  class Collapse;
  namespace details
  { 
    template< class T_FUNCTION,
              bool arity = boost::function_traits<T_FUNCTION> :: arity >
       class FreeFunction;
  }
  //! \endcond

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
      Base() : basis() { coefs.resize( basis.size() ); }
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

    protected:
      //! Links basis functions.
      t_GroupOp groupop;
      //! Links scalars to basis functions.
      t_ScalarOp scalarop;

    protected:
      //! For alternating least-square purposes.
      template< class T_ITIN, class T_ITOUT>
        t_Return expand( T_ITIN _first, T_ITIN _last, T_ITOUT _out ) const;
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


  /** \brief Collapses a sum of separable function into Fitting::Allsq format.
   *  \details This flavor keeps track of computed basis functions
   *           (\f$g_{d,i}^{(r)}\f$).  This is quite a complex process. A sum
   *           of separable functions is expressed as follows: 
   *    \f[
   *        F(\mathbf{x}) = \sum_r \prod_d \sum_i \lambda_{d,i}^{(r)}
   *                        g_{i,n}^{(r)}(x_i),
   *    \f]
   *    We will keep track of the r, i, and n indices. Furthermore, the least
   *    square-fit method fits to  \e o observables. It expects a "2d" input matrix
   *    I[d, (r,i) ], where the slowest running index is the dimension \e d.
   *    I[d] should be a vecrtor container type of scalar values. These scalar
   *    values are orderer with  \e i the fastest running index and \e r the
   *    slowest. The method also expects a functor which can create the matrix
   *    \e A (as in \e Ax = \e b ) for each specific dimension \e d. The type
   *    of A is a vector of scalars. The fastest running index is \e i,
   *    followed by \e r, and finally \e o. 
   *    \tparam T_FUNCTION must be a Function< Summand< T_BASIS > > 
   **/
  template< class T_FUNCTION >
  class Collapse
  {
    public:
      //! Type of the sum of separable functions to collapse.
      typedef T_FUNCTION t_Function;
      //! Wether to update the coefficients between each dimension or not.
      bool do_update;
      //! Reference to the function to fit.
      T_FUNCTION &function;

      //! Constructor
      Collapse   ( t_Function &_function )
                  : do_update(false), is_initialized(false), D(0), nb_obs(0),
                    nb_ranks(0), function( _function ) {}

      //! \brief Constructs the matrix \a A for dimension \a _dim. 
      //! \details Note that depending \a _coef are used only for dim == 0,
      //!          and/or do_update == true.
      //! \tparam T_MATRIX is a vector, e.~g. a flattened matrix. This class is
      //!                  meant for %GSL which expects this type of memory
      //!                  setup for matrices. 
      //! \tparam T_VECTORS is a vector of vectors or similar. 
      //! \param[out] _A will contain the collapsed matrix for the linear
      //!                least-square fit. \a _A[ (o,r,i ) ].
      //! \param[in] The dimension for which to create \a _A.
      //! \param[in] _coefs contains the coefficients for \e all dimensions.
      //!                   \a _coefs[d, (r,i) ]
      template< class T_MATRIX, class T_VECTORS >
      void operator()( T_MATRIX &_A, types::t_unsigned _dim, const T_VECTORS &_coefs );

      //! \brief Constructs the completely expanded matrix.
      //! \tparams T_VECTORS is a vector of vectors.
      //! \params[in] _x contains the input to sum of separable function in the
      //!                format \a _x[ o,d ].
      template< class T_VECTORS > void init( const T_VECTORS &_x );
      //! Resets collapse functor. Clears memory.
      void reset();
      //! Creates a collection of random coefficients.
      //! \details A vector of vectors of correct dimensions is created but not
      //!          necessarily initialized.
      //! \tparams should be a vector of vectors.
      //! \params[out] _coefs creates \a _coefs[d, (r,i) ] uninitialized.
      template< class T_VECTORS > void create_coefs( T_VECTORS &_coefs ) const;
      //! Assigns solution coefficients to Collapse::function.
      //! \tparam T_VECTORS is a vector of vectors or similar. 
      //! \param[in] _solution contains the coefficients for \e all dimensions.
      //!                      \a _solution[d, (r,i) ]
      template< class T_VECTORS > void reassign( const T_VECTORS &_solution ) const;

      __DODEBUGCODE
      ( 
        //! \cond
        template<class T_VECTORS> void print_funcs( const T_VECTORS &_x ) const; 
        //! \endcond
      )

    protected:
      //! Initializes Collapse::factors.
      //! \tparam T_VECTORS is a vector of vectors or similar. 
      //! \param[in] _coefs contains the coefficients for \e all dimensions.
      //!                   \a _coefs[d, (r,i) ]
      template< class T_VECTORS >
        void initialize_factors( const T_VECTORS &_coefs );
      //! Updates Collapse::factors.
      //! \tparam T_VECTORS is a vector of vectors or similar. 
      //! \param[in] The dimension for which to create \a _A. Passed from
      //!            Collapse::operator()().
      //! \param[in] _coefs contains the coefficients for \e all dimensions.
      //!                   \a _coefs[d, (r,i) ].
      template< class T_VECTORS >
        void update_factors( types::t_unsigned _dim, const T_VECTORS &_coefs );

      //! False if unitialized.
      bool is_initialized;
      //! Maximum dimension.
      types::t_unsigned D;
      //! Number of observables.
      types::t_unsigned nb_obs;
      //! Number of ranks.
      types::t_unsigned nb_ranks;
      //! Return type of the function
      typedef typename T_FUNCTION :: t_Return t_Type;
      //! \brief Type of the matrix containing expanded function elements.
      typedef std::vector< std::vector< t_Type > > t_Expanded;
      //! A matrix with all expanded function elements.
      //! \details expanded[ (o,d), (r,i) ]. The type is a vector of vectors.
      //!          The fastest-running \e internal index is \e i. The
      //!          fastest-running \e external index is d.
      //!          \f$ \text{expanded}[(o,d), (r,i)] = g_{d,i}^{(r)}(x_d^{(o)})\f$.
      t_Expanded expanded;
      //! Type of the factors.
      typedef std::vector< std::vector< t_Type > > t_Factors;
      /** \brief A matrix wich contains the factors of the separable functions.
       *  \details factors[(o,r), d]. A vector of vectors. The \e internal
       *           index is \e d. The fastest-running \e external
       *           index is \e r.
       *           \f$
       *               \text{factors}[d, (r,i)] = \sum_i
       *               \lambda_{d,i}^{(r)}g_{d,i}^{(r)}(x_d^{(o)})
       *           \f$. **/
      t_Factors factors;
      //! Type of the sizes.
      typedef std::vector< std::vector< types::t_unsigned > > t_Sizes;
      //! \brief Sizes of the basis per dimension and rank.
      //! \details sizes[ d, r ] = max i \@ (d,r). r is the fastest running
      //!          index.
      t_Sizes sizes;
  };

} // end of Separable namespace

#include "separable.impl.h"

#endif
