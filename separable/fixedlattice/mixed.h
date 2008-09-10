//
//  Version: $Id$
//
#ifndef _CE_MIXED_H_
#define _CE_MIXED_H_

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

namespace Traits 
{
  namespace CE
  {
    //! \brief Traits for a combination of Cluster Expansion with Sum of
    //!       Separable Functions.
    template< class T_COLLAPSE, class T_CEBASE >
     struct MixedApproach
     {
       protected:
         //! Original separables.
         typedef typename T_COLTRAITS :: t_Separables t_OrigSeparables;

       public:
         //! Original separable traits.
         typedef typename t_OrigSeparables t_Traits t_OrigSepTraits;
         //! Type of the separable function.
         typedef typename SeparablesWithMatrixRange<t_OrigSeparables>
                                                   :: other t_Separables;
         //! New Collapse functor.
         typedef typename CollapseWithNewSeparables<t_Separables> :: other t_Collapse;
         //! Type of the configuration matrix.
         typedef typename t_Collapse :: t_Configurations t_Configurations;
         //! Type of the Mapping.
         typedef typename t_Collapse :: t_Mapping t_Mapping;
         //! The Regularization policy.
         typedef typename t_Collapse :: t_RegPolicy t_RegPolicy;
         //! Type of the Update Policy.
         typedef typename t_Collapse :: t_UpdatePolicy t_UpdatePolicy;
         //! CE fit base.
         typedef T_CEBASE t_CEBase;
     };

  } // end of CE namespace 

} // end of Traits namespace
namespace CE
{

  //! \brief Combination of Cluster Expansion with Sum of Separable Functions.
  template< class T_TRAITS >
    class MixedApproach
    {
      template< class TT_TRAITS > friend class MixedApproach;
      public:
        //! Type of the traits.
        typedef T_TRAITS t_Traits;
        //! Type of the separable collapse functor.
        typedef typename t_Traits :: t_Collapse t_Collapse;
        //! Type of the CE fitting base class.
        typedef typename t_Traits :: t_CEBase  t_CEFit;
        //! Type of the separable function.
        typedef typename t_Traits :: t_Separables  t_Separables;
        //! Type of the configuration matricex.
        typedef typename t_Traits :: t_Configurations t_Configurations;
        //! \brief Type of the matrix range.
        //! \details Necessary interface for minimizer.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Matrix t_Matrix;
        //! Type of the vectors.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Vector t_Vector;
        //! Type of the classes of equivalent clusters.
        typedef std::vector< std::vector< ::CE::Cluster > > t_Clusters;


        //! Constructor.
        MixedApproach() 
          { clusters_.reset( new t_Clusters ); coefficients_.reset( new t_Matrix ); }
        //! Constructor.
        template< class TT_TRAITS >
          MixedApproach   ( const MixedApproach< TT_TRAITS > &_c )
                        : separables_( _c.separables_ ), 
                          clusters_( _c.clusters_ ), 
                          coefficients_( _c.coefficients_ ),
                          collapse_( _c.collapse_ ),
                          cefit_( _c.cefit_ )
            { init( _c.separables().ranks(), _c.separables().dimensions() ); }
        //! Destructor.
        ~MixedApproach() {}

        //! Returns a reference to the CE fitting base.
        t_CEFit& cefit() { return cefit_; }
        //! Returns a constant reference to the CE fitting base.
        const t_CEFit& cefit() const { return cefit_; }
        //! Returns a reference to the CE fitting base.
        t_Collapse& collapse() { return collapse_; }
        //! Returns a constant reference to the CE fitting base.
        const t_Collapse& collapse() const { return collapse_; }


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
        void reset() { collapse().reset(); }

        //! Reference to configuration matrix.
        const size_t nbconfs() const { return collapse().nbconfs(); }
        
        //! Returns the number of dimensions.
        size_t dimensions() const { return collapse().dimensions(); }
        //! Returns the number of degrees of liberty (per dimension).
        size_t dof() const { return collapse().dof() + cefit().dof(); }
       
        //! Initializes the mixed approach.
        void init( size_t _ranks, size_t _dims );

        //! Returns a reference to the coefficients.
        t_Matrix& coefficients()
        {
          __ASSERT( not coefficients_.get(), "Empty smart pointer.\n" ) 
          return *coefficients_; 
        }
        //! Returns a constant reference to the coefficients.
        const t_Matrix& coefficients() const 
        {
          __ASSERT( not coefficients_.get(), "Empty smart pointer.\n" ) 
          return *coefficients_; 
        }
        //! Returns a reference to the mapping.
        typename t_Collapse::t_Mapping& mapping() { return collapse().mapping(); }
        //! Returns a reference to the mapping.
        const typename t_Collapse::t_Mapping& mapping() const
          { return collapse().mapping(); }
        //! Returns a reference to the separable function;
        t_Separables& separables() { return separables_; }
        //! Returns a constant reference to the separable function;
        const t_Separables& separables() const { return separables_; }
        //! Returns a reference to the clusters.
        t_Clusters& clusters()
        {
          __ASSERT( not clusters_.get(), "Empty smart pointer.\n" ) 
          return *clusters_; 
        }
        //! Returns a constant reference to the clusters.
        const t_Clusters& clusters() const 
        {
          __ASSERT( not clusters_.get(), "Empty smart pointer.\n" ) 
          return *clusters_; 
        }
        //! Reassigns clusters.
        void reassign();
        //! Randomizes both cluster energies and ecis.
        void randomize( typename t_Vector :: value_type _howrandom );

      protected:
        //! Creates the _A and _b matrices for fitting.
        template< class T_MATRIX, class T_VECTOR >
          void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );

        //! The separable function.
        t_Separables separables_;

        //! The coefficients.
        boost::shared_ptr<t_Matrix> coefficients_;

        //! The classes of equivalent clusters for mixing.
        boost::shared_ptr<t_Clusters> clusters_;

        //! The collapse functor.
        t_Collapse collapse_;

        //! The CE fitting functor.
        t_CEFit cefit_;

        //! Current dimension we are working on.
        size_t dim;
    };

  //! Prints mixed-approach description to a stream.
  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream,
                            const MixedApproach<T_TRAITS> &_col );

  //! \brief Removes cluster classes for wich at least one custer is contained
  //!        within \a _positions and the origin.
  template< class T_POSITIONS, class T_CLUSTERS >
    void remove_contained_clusters( const T_POSITIONS &_positions, T_CLUSTERS &_clusters );
} // end of CE namspace

#include "mixed.impl.h"

#endif
