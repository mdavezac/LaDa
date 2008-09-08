//
//  Version: $Id$
//
#ifndef _CE_MANY_H_
#define _CE_MANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include<iostream>
#include<vector>
#include<list>

#include <opt/types.h>
#include <opt/debug.h>

namespace CE
{
  template< class T_TRAITS >
    class Many 
    {
      public:
        //! Type of the traits.
        typedef T_TRAITS t_Traits;
        //! Type of the separable collapse functor.
        typedef typename t_Traits :: t_Collapse t_Collapse;
        //! \brief Type of the matrix range.
        //! \details Necessary interface for minimizer.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Matrix t_Matrix;
        //! Type of the vectors.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Vector t_Vector;
        //! Type of the container of separables.
        boost::ptr_list< t_Collapse > t_Collapses;
        //! Type of the container of separables.
        boost::ptr_list< typename t_Traits :: t_Separables > t_CtnrSeparables;


        //! Constructor.
        Many() {}
        //! Destructor.
        ~Many() {}

        //! Creates the fitting matrix and target vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &_A, T_VECTOR &_b,
                           types::t_unsigned _dim );
        //! Evaluates square errors.
        opt::ErrorTuple evaluate() const;
        //! Predicts target value of a structure.
        typename t_Matrix :: value_type evaluate( size_t _n ) const;

        //! Updates the separable and copies the eci from column 0 to all other columns.
        void update_all();
        //! Updates the separable and copies the eci from column d to column 0.
        void update( types::t_unsigned _d );
        //! Resets collapse functor.
        void reset() { t_ColBase::reset(); }

        //! Returns the number of dimensions.
        size_t dimensions() const;
        //! Returns the number of degrees of liberty (per dimension).
        size_t dof() const;
       
        //! Initializes the mixed approach.
        void init( size_t _ranks, size_t _dims );

        //! Randomizes both cluster energies and ecis.
        void randomize( typename t_Vector :: value_type _howrandom );

      protected:
        //! Returns the number of degrees of liberty for current dimension.
        size_t current_dof() const;
        //! Creates the _A and _b matrices for fitting.
        template< class T_MATRIX, class T_VECTOR >
          void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );
        //! The container of separable functions.
        t_CtnrSeparables separables_;
        //! The collapse functor associated with the separable functions.
        t_Collapses collapses_;
        //! Current dimension being updated.
        size_t dim;
    };

  //! Prints mixed-approach description to a stream.
  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream,
                            const Many<T_TRAITS> &_col );

} // end of CE namespace

#include "many.impl.h"

#endif
