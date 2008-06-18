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
  /** \brief A sum of separable functions.
   * \details The separable function \f$F(\mathbf{x})\f$ acting upon vector
   *          \f$\mathbf{x}\f$ and returning a scalar can be defined as 
   *    \f[
   *        F(\mathbf{x}) = \sum_r \prod_i \sum_n \lambda_{i,n}^{(r)}
   *                        g_{i,n}^{(r)}(x_i),
   *    \f]
   *          with the sum over \e r running over the ranks of the separable
   *          functions, the product over \e i are the separable functions
   *          proper, and the sum_n is an expansion of the factors ver some
   *          family of 1d-functions.
   * \tparam T_BASIS is a container of 1d functions. These functions should
   *         zero order evaluation via a functor call, and grdient evaluation
   *         via a t_Return gradient( t_Arg ) member function. **/
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

  template< class T_FUNCTION >
  class Collapse
  {
    public:
      typedef T_FUNCTION t_Function;

      //! Constructor
      Collapse   ( T_FUNCTION &_function )
               : is_initialized(false), D(0), nb_obs(0),
                 nb_ranks(0), function( _function ) {}

      template< class T_MATRIX, class T_VECTORS >
      void operator()( T_MATRIX &_A, types::t_unsigned _dim, T_VECTORS &_coefs )

      //! Constructs the completely expanded matrix.
      void init();

    protected:
      //! False if unitialized.
      bool is_initialized;
      //! Maximum dimension.
      types::t_unsigned D;
      //! Number of observables.
      types::t_unsigned nb_obs;
      //! Number of ranks.
      types::t_unsigned nb_ranks;
      //! Reference to the function to fit.
      T_FUNCTION &function;
      //! Return type of the function
      typedef typename T_FUNCTION :: t_Return t_Type;
      //! \brief Type of the matrix containing expanded function elements.
      //! \details This is, more or less, a 2d matrix. Each row contains the
      //!          \f$ g_{i,n}^{(r)}(x_i^{(u)}) \f$ described in Separable::Function for a
      //!          specific dimension \e i and an observable
      //!          \f$\mathbs{x}^{(u)}\f$. The rows encoded in \e n major and
      //!          \e r minor. The columns are encoded in \e i major, \e u minor.
      typedef std::vector< std::vector< t_Type > > t_Expanded;
      //! A matrix with all expanded function elements.
      t_Expanded expanded;
      //! Type of the factors.
      typedef std::vector< std::vector< t_Type > > t_Factors;
      /** \brief A matrix wich is exanded only up to \e i.
       *  \details The sum over \e n in the expression of Separable::Function
       *           is not expanded. Otherwise, the encoding is similar to the
       *           one described in Separable :: Collapse :: t_Expanded, with
       *           \f$ \phi_i^{(r)}( x_i^{(u)} ) = \sum_n
       *           \lambda_{i,n}^{(r)}g_{i,n}^{(r)}( x_i^{(u)} )\f$. Note that
       *           this time the coefficients are contained in the matrix.
       **/
      t_Factors factors;
      //! Type of the sizes.
      typedef std::vector< std::vector< types::t_unsigned > > t_Sizes;
      //! Sizes of the basis per dimension and rank.
      t_Sizes sizes;
  }

  template<class T_FUNCTION> template< class T_VECTORS >
  void Collapse<T_FUNCTION> :: init( const T_VECTORS &_x ) 
  {
    // finds maximum size of separable function.
    {
      using namespace boost::lambda;
      typename t_Function::t_Basis::const_iterator i_max;
      i_max = std::max_element( _x.basis.begin(), x.basis.end(), 
                                bind( &T_FUNCTION::t_Basis::value_type::size, _1 ) 
                                < 
                                bind( &T_FUNCTION::t_Basis::value_type::size, _2 )  );
      D = i_max->size();
    }
                      
    nb_obs = _x.size();
    nb_ranks = function.basis.size();
    sizes.resize( D );
    expanded.reserve( nb_obs * D );
    typename T_VECTORS::const_iterator i_x = _x.begin(); 
    typename T_VECTORS::const_iterator i_x_end = _x.end(); 
    for(; i_x != i_x_end; ++i_x ) // loop over observables.
    {
      typename T_VECTORS::value_type::const_iterator i_var = i_x->begin();
      for( size_t dim = 0; dim < D; ++dim, ++i_var )
      {
        __ASSERT( i_var != i_x->end(), "Inconsistent sizes.\n" )
        __ASSERT( nb_ranks != function.basis[dim].size(),
                  "Inconsistent number of ranks.\n" )

        expanded.resize( expanded.size() + 1 );

        typename t_Expanded :: value_type &row = expanded.back();
        typedef typename t_Function::t_Basis::const_iterator t_const_iterator;
        t_const_iterator i_rank = function.basis[dim].begin();
        t_const_iterator i_rank_end = function.basis[dim].end();
        for(; i_rank != i_rank_end; ++i_rank )
        {
          if( i_rank->size() < dim ) continue;
          sizes[dim].push_back( i_rank->basis.size() );
          typedef typename t_Function::t_Basis::value_type
                                     ::t_Basis::const_iterator t_citerator;
          t_citerator i_func = i_rank->basis.begin();
          t_citerator i_func_end = i_rank->basis.end();
          for(; i_func != i_func_end; ++i_func )
            row.push_back( (*i_func)(*i_var) );
        }
      }
    }
  }

  template< class T_MATRIX, class T_VECTORS >
  void Collapse<T_FUNCTION>::operator()( T_MATRIX &_A, types::t_unsigned _dim,
                                         const T_VECTORS &_coefs )
  {
    if( ( not is_initialized ) or ( ( not do_update ) and _dim == 0 ) )
      initialize_factors( _coefs );
    else if ( do_update ) update_factors( _dim, _coefs );
    is_initialized = true;   

    _A.resize( nb_obs );
    typename T_MATRIX :: iterator i_A = _A.begin();
    typename T_MATRIX :: iterator i_A_end = _A.end();
    typename t_Expanded :: const_iterator i_ex = expanded.begin() + dim;
    for(; i_A != i_A_end; i_ex += dim ) // loops over observables
    {
      t_Size :: const_iterator i_size = sizes[dim].begin();
      t_Size :: const_iterator i_size_end = sizes[dim].end();
      t_Factors :: const_iterator i_facs = factors.begin();
      __ASSERT( factors.size() != sizes[dim].size(),
                "Inconsistent sizes.\n" )
      typename t_Expanded :: value_type const_iterator i_var = i_ex->begin();
      for(; i_size != i_size_end; ++i_size, ++i_facs ) // loops over ranks.
      {
        // Computes factor for this rank.
        typename t_Factors :: value_type :: const_iterator i_fac = i_facs.begin();
        typename t_Factors :: value_type :: const_iterator i_fac_end = i_facs.end();
        typename t_Factors ::value_type :: value_type factor(1);
        for(types::t_unsigned i = 0; i_fac != i_fac_end; ++i_fac, ++i )
          if( i != Dim ) factor *= (*i_fac);

        // loops over basis set of factor function of current dimension and rank.
        for( types::t_unsigned i = *i_size; i > 0; --i, ++i_var, ++i_A )
        {
          __ASSERT( i_var != expanded.end(), "Inconsistent sizes.\n" )
          __ASSERT( i_A != _A.end(), "Inconsistent sizes.\n" )
          *i_A = (*i_var) * factor;
        }
      }
    }
  }

  template< class T_MATRIX, class T_VECTORS >
  void Collapse<T_FUNCTION>::initialize_factors( const T_VECTORS &_coefs )
  {
    __ASSERT( nb_ranks != function.basis.size(),
              "Inconsistent rank size.\n" )
    factors.resize( nb_ranks );
    typename t_Function :: t_Basis :: const_iterator i_rank = function.basis.begin();
    typename t_Factors :: iterator i_facs = factors.begin();
    typename t_Factors :: iterator i_facs_end = factors.end();
    typename t_Expanded :: const_iterator i_ex = expanded.begin();
    for(; i_facs != i_facs_end; ++i_rank, ++i_facs )
    {
      i_facs->resize( i_rank->basis.size() );
      typename t_Factors :: value_type :: iterator i_fac = i_fac->begin();
      typename t_Factors :: value_type :: iterator i_fac_end = i_fac->end();
      for(types::t_unsigned dim = 0; i_fac != i_fac_end; ++i_fac, ++dim )
      {
        *i_fac = t_Type(1);
      }
    }
  }
}

#include "separable.impl.h"

#endif
