//
//  Version: $Id$
//
#ifndef _CE_MANY_H_
#define _CE_MANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/ptr_container/ptr_list.hpp>
#include<boost/shared_ptr.hpp>
#include<boost/mpl/apply.hpp>
#include<boost/mpl/transform.hpp>
#include<boost/fusion/adapted/mpl.hpp>
#include<boost/fusion/include/mpl.hpp>
#include<boost/fusion/container/list.hpp>


#include<iostream>
#include<vector>
#include<list>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/indirection.h>

//! \cond 
namespace CE { template< class T_TRAITS > class Many; }
//! \endcond

namespace Traits
{
  namespace CE
  {
    //! Traits of a "Many" collapse functor.
    template< class T_LIST_OF_SEPARABLES,
              class T_LIST_OF_COLLAPSES,
              class T_MAPPING = ::CE::Mapping::Basic,
              class T_COEFFICIENTS = boost::numeric::ublas::matrix<types::t_real>,
              class T_VECTOR
                = boost::numeric::ublas::vector<typename T_COEFFICIENTS::value_type> >
    struct Many 
    {
      protected:
        //! Transforms a separable into 
        typedef SeparablesWithMatrixRange< boost::mpl::_ > make_sep;
        typedef CollapseWithNewSeparables< boost::mpl::_ > make_col;
        //! Transforms an type into the type of a container of vectors.
        typedef boost::ptr_list< boost::mpl::_ > :: other make_ptrlist;
        //! List of types in a mpl::list type.
        typedef boost::mpl::transform
                < 
                  T_LIST_OF_SEPARABLES,
                  boost::mpl::apply
                  < 
                    make_ptrlist,
                    make_sep,
                  > :: type 
                >  t_MPLListOfSeparables;
        //! List of types in a mpl::list type.
        typedef mpl::transform< T_LIST_OF_COLLAPSES, make_ptrlist > :: type 
                t_MPLListOfCollapses;
      public:
        //! Tuple of containers of separables.
        typedef t_MPLListOfSeparables t_ListsOfSeparables;
        //! Tuple of containers of collapses.
        typedef t_MPLListOfSeparables t_ListsOfCollapses;
        //! Type of the Mapping.
        typedef T_MAPPING t_Mapping;
        //! Type of the coefficients.
        typedef T_COEFFICIENTS t_Coefficients;
        //! Type of the vectors.
        typedef T_VECTOR t_Vector;
    };
  }

} // end of traits namespace.

namespace CE
{
#define VOID_MFUNC( FunctionName ) Void_Unary_Apply ## FunctionName
#define VOID_MFUNC_DECLARE( FunctionName ) \
  struct ACC_METAFUNCTIONNAME( FunctionName )\
  {\
    template< class T > void operator()( T &_t ) \
    {\
      foreach( const typename T::value_type &_val, _t )\
        _val.FunctionName(); \
      return _r;\
    }\
  }; 
#define ACC_MFUNC( FunctionName ) AccApply ## FunctionName
#define ACC_MFUNC_DECLARE( FunctionName ) \
  struct ACC_METAFUNCTIONNAME( FunctionName )\
  {\
    template< class T > size_t operator()( T &_t, size_t _r )\
    {\
      foreach( const typename T::value_type &_val, _t )\
        _r += _val.FunctionName(); \
      return _r;\
    }\
  }; 
#define U_MFUNC( FunctionName, Arg ) UApply ## FunctionName( Arg )
#define U_MFUNC_DECLARE( FunctionName, Type ) \
  struct U_Apply ## FunctionName \
  {\
    Type &arg;\
    U_Apply ## FunctionName ( const Type &_t ) : arg( _t ) {}\
    U_Apply ## FunctionName ( const U_Apply ## Functioname &_c ) : arg( _c.arg ) {}\
    template< class T > void operator()( T &_t )\
    {\
      foreach( const typename T::value_type &_val, _t )\
        _val.FunctionName( arg );\
    }\
  };

  namespace details
  {
    //! Redefine separables traits such that coefficients are a range.
    template< class T_TRAITS, class T_COEFFICIENTS >
    struct SepTraits
    {
      //! Mapping traits.
      typedef typename T_TRAITS :: t_Mapping t_Mapping;
      //! Policy traits.
      typedef typename T_TRAITS :: t_Policy t_Policy;
      //! Range coefficients traits.
      typedef ::CE::Policy::MatrixRangeCoefficients t_Coefficients;
      //! Range coefficients traits.
      typedef typename t_Coefficients :: t_Matrix t_Matrix;
      //! Range coefficients traits.
      typedef typename T_TRAITS :: t_Vector t_Vector;
    };
    //! Redefines collapse traits with new separables.
    template< class T_SEPARABLES, class T_TRAITS >
    struct ColTraits
    {
      //! Type of the configuration matrix.
      typedef typename T_TRAITS :: t_Configurations t_Configurations;
      //! Type of the Mapping.
      typedef T_SEPARABLES t_Separables;
      //! Type of the Mapping.
      typedef typename T_TRAITS :: t_Mapping t_Mapping;
      //! Type of the Regulations Policy
      typedef typename T_TRAITS :: t_RegPolicy t_RegPolicy;
      //! Type of the Policy.
      typedef typename T_TRAITS :: t_UpdatePolicy t_UpdatePolicy;
    };
  }

  template< class T_TRAITS >
    class Many 
    {
      public:
        //! Type of the traits.
        typedef T_TRAITS t_Traits;
        //! Type of the of separables function.
        typedef opt::IndirectionBase  t_Separables;
        //! Type of the of separables function.
        typedef opt::IndirectionBase  t_Collapse;
        //! \brief Type of the matrix coefficients.
        typedef typename t_Traits :: t_Coefficients t_Coefficients;
        //! \brief Type of the matrix range.
        //! \details Necessary interface for minimizer.
        typedef typename t_Traits :: t_Coefficients t_Matrix;
        //! Type of the vectors.
        typedef typename t_Traits :: t_Vector t_Vector;
        //! Type of the container of separables.
        typedef boost::shared_ptr< t_ListsOfCollapse > t_ListsOfCollapses;
        //! Type of the container of separables.
        typedef boost::shared_ptr< t_ListsOfSeparables > t_ListsOfSeparables;
        //! Type of the general mapping.
        typedef typename t_Traits :: t_Mapping t_Mapping;


        //! Constructor.
        Many() : separables_( new t_CtnrSeparables ),
                 collapses_( new t_Collapses ), dim(0) {}
        //! Copy Constructor.
        Many( const Many& _c ) : separables_( _c.separables_ ),
                                 collapses_( _c.collapses ),
                                 dim( _c.dim ), mapping_( _c.mapping_ ),
                                 coefficients_( _c.coefficients_ ) {}
        //! Destructor.
        ~Many() {}

        //! Creates the fitting matrix and target vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &_A, T_VECTOR &_b,
                           types::t_unsigned _dim );
        //! Evaluates square errors.
        opt::ErrorTuple evaluate() const;
        //! Predicts target value of a structure.
        typename t_Matrix :: value_type evaluate( size_t _n ) const
         { return boost::fusion::accumulate( *collapses_, ApplyEvaluateOne(_n) ); }


        //! \brief Updates the separable and copies the eci from column 0 to all
        //!        other columns.
        void update_all();
         { return boost::fusion::for_each( *collapses_, VOID_MFUNC(update_all)() ); }
        //! Updates the separable and copies the eci from column d to column 0.
        void update( types::t_unsigned _d )
         { return boost::fusion::for_each( *collapses_, U_MFUNC(update)(_d) ); }
        //! Resets collapse functor.
        void reset();
         { return boost::fusion::accumulate( *collapses_, VOID_MFUNC(reset)() ); }

        //! Returns the number of dimensions.
        size_t dimensions() const
         { return boost::fusion::accumulate( *collapses_, ACC_MFUNC(dimensions)() ); }
        //! Returns the number of degrees of liberty (per dimension).
        size_t dof() const
         { return boost::fusion::accumulate( *collapses_, ACC_MFUNC(dof)() ); }
        //! Returns the number of configurations.
        size_t nbconfs() const
         { return boost::fusion::accumulate( *collapses_, ACC_MFUNC(nbconfs)() ); }
       
        //! Randomizes both cluster energies and ecis.
        void randomize( typename t_Vector :: value_type _howrandom )
         { return boost::fusion::for_each( *collapses_, U_MFUNC(randomize)(_howrandom) ); }

        //! Add new collapse and separables.
        template< class T_COLLAPSE, class T_SEPARABLES > size_t add_as_is();
        //! Add new collapse and separables.
        template< class T_COLLAPSE, class T_SEPARABLES > size_t wrap_n_add();
        
        //! Returns reference to nth separable function.
        t_Separables* separables( size_t _n ) { return (*separables_)[_n].self(); }
        //! Returns constant reference to nth separable function.
        const t_Separables* separables( size_t _n ) const
          { return (*separables_)[_n].self(); }
        //! Returns reference to nth collapse functor.
        t_Collapse* collapse( size_t _n ) { return (*collapses_)[_n].self(); }
        //! Returns constant reference to nth collapse functor.
        t_Collapse* const collapse( size_t _n ) const
          { return (*collapses_)[_n].self(); }
        //! Returns the number of collapse and separables functions.
        size_t size() const { return collapses_->size(); }
        //! Returns a reference to the mapping.
        t_Mapping mapping() { return mapping_; }
        //! Returns a constant reference to the mapping.
        const t_Mapping mapping() const { return mapping_; }
        //! Initializes a collapse with rank and size.
        void init( size_t _index, size_t _rank, size_t _dimensions );

      protected:
        //! Returns the number of degrees of liberty for current dimension.
        size_t current_dof() const
         { return boost::fusion::accumulate( *collapses_, ApplyCurrentDof(dim) ); }
        //! Creates the _A and _b matrices for fitting.
        template< class T_MATRIX, class T_VECTOR >
          void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );

        //! The container of separable functions.
        t_ListsOfSeparables separables_;
        //! The collapse functor associated with the separable functions.
        t_ListsOfCollapses collapses_;
        //! Current dimension being updated.
        size_t dim;
        //! The mapping to the structures ( e.g. leave-one-out, leave-many-out )
        t_Mapping mapping_;
        //! The coefficienst.
        t_Coefficients coefficients_;

        // metafunctions.
        //! \cond
        ACC_MFUNC_DECLARE(dimensions)
        ACC_MFUNC_DECLARE(dof)
        ACC_MFUNC_DECLARE(nbconfs)
        struct ApplyCurrentDof;
        VOID_MFUNC_DECLARE(reset);
        template< class T_VECTOR > struct ApplyCreateAnB;
        template< class T_MATRIX, class T_VECTOR > struct ApplyRegularization;
        U_MFUNC_DECLARE(randomize, typename t_Vector::value_type );
        VOID_MFUNC_DECLARE(update_all);
        U_MFUNC_DECLARE(update, types::t_unsigned _d);
        struct ApplyEvaluateOne;
        //! \endcond
    };

  //! Prints mixed-approach description to a stream.
  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream, const Many<T_TRAITS> &_col );

  //! Initializes a Many separable function depending on string input and structures.
  template< class T_STRUCTURES, class T_TRAITS >
   void init_many_collapses( const std::string &_desc, size_t _rank,
                             types::t_real _lambda, const T_STRUCTURES &_structures,
                             Many<T_TRAITS> &_many );

} // end of CE namespace

#include "many.impl.h"

#define VOID_MFUNC
#define VOID_MFUNC_DECLARE
#define ACC_MFUNC
#define ACC_MFUNC_DECLARE
#define U_MFUNC
#define U_MFUNC_DECLARE

#endif
