#ifndef _CE_MIXED_H_
#define _CE_MIXED_H_

#include "LaDaConfig.h"

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include<iostream>
#include<vector>
#include<list>

#include <misc/types.h>
#include <opt/debug.h>

namespace LaDa
{
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
           typedef typename T_COLLAPSE :: t_Traits :: t_Separables t_OrigSeparables;

         public:
           //! Original Collapse type.
           typedef T_COLLAPSE t_OrigCollapse;
           //! Original separable traits.
           typedef typename t_OrigSeparables :: t_Traits t_OrigSepTraits;
           //! Type of the separable function.
           typedef typename SeparablesWithMatrixRange<t_OrigSeparables>
                                                     :: type t_Separables;
           //! New Collapse functor.
           typedef typename CollapseWithNewSeparables<t_OrigCollapse, t_Separables>
                                :: type t_Collapse;
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

           //! Rebinds traits.
           template< class TT_COLLAPSE, class TT_CEBASE >
           struct rebind
           {
             //! Result type.
             typedef MixedApproach< TT_COLLAPSE, TT_CEBASE >  type;
           };

  //       //! Type of the coefficients range.
  //       typedef typename t_Traits :: t_OrigSepTraits :: t_Matrix t_Coefficients;
       };

      //! Rebinds collapse functor with new mapping.
      template< class T_MIXED, class T_MAPPING >
      struct MixedApproachWithNewMapping
      {
        protected:
          //! Original type.
          typedef T_MIXED t_Mixed;
          //! Original traits.
          typedef typename t_Mixed :: t_Traits t_Traits;
          //! Original CE type`.
          typedef typename t_Traits :: t_CEBase t_CEBase;
          //! New Collapse. 
          typedef typename t_Traits :: t_OrigCollapse t_Collapse;
          //! Original Collapse traits.
          typedef typename t_Collapse :: t_Traits t_CollapseTraits;
          //! New Collapse Traits
          typedef typename t_CollapseTraits ::template rebind
                  <
                    typename t_CollapseTraits :: t_Separables, 
                    T_MAPPING,
                    typename t_CollapseTraits :: t_RegPolicy, 
                    typename t_CollapseTraits :: t_Configurations,
                    typename t_CollapseTraits :: t_UpdatePolicy
                  > :: type t_NewCollapseTraits;
          //! New Collapse. 
          typedef typename t_Collapse ::template rebind
                  <
                    t_NewCollapseTraits
                  > :: type t_NewCollapse;
          //! New traits
          typedef typename t_Traits ::template rebind
                  < 
                    t_NewCollapse, 
                    t_CEBase
                  > :: type t_NewTraits;
        public:
          //! Resulting type.
          typedef typename t_Mixed ::template rebind
                  <
                    t_NewTraits 
                  > :: type type;
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
          //! Rebinds MixedApproach.
          template< class TT_TRAITS >
          struct rebind
          {
            //! Result type.
            typedef MixedApproach< TT_TRAITS >  type;
          };
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
          typedef std::vector< std::vector< ::LaDa::CE::Cluster > > t_Clusters;


          //! Constructor.
          MixedApproach() {}
          //! Constructor.
          MixedApproach   ( const MixedApproach &_c )
                        : clusters_( _c.clusters_ ), 
                          coefficients_( _c.coefficients_ ),
                          collapse_( _c.collapse_ ),
                          cefit_( _c.cefit_ )
            { init( _c.separables().ranks(), _c.separables().dimensions() ); }
          //! Destructor.
          ~MixedApproach() {}
          //! Assignement operator.
          template< class TT_TRAITS >
            void operator=( const MixedApproach<TT_TRAITS> &_c );

          //! Returns a reference to the CE fitting base.
          t_CEFit& cefit() { return cefit_(); }
          //! Returns a constant reference to the CE fitting base.
          const t_CEFit& cefit() const { return cefit_(); }
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
          t_Matrix& coefficients() { return coefficients_(); }
          //! Returns a constant reference to the coefficients.
          const t_Matrix& coefficients() const { return coefficients_(); }
          //! Returns a reference to the mapping.
          typename t_Collapse::t_Mapping& mapping() { return collapse().mapping(); }
          //! Returns a reference to the mapping.
          const typename t_Collapse::t_Mapping& mapping() const
            { return collapse().mapping(); }
          //! Returns a reference to the separable function;
          t_Separables& separables() { return collapse().separables(); }
          //! Returns a constant reference to the separable function;
          const t_Separables& separables() const { return collapse().separables(); }
          //! Returns a reference to the clusters.
          t_Clusters& clusters() { return clusters_(); }
          //! Returns a constant reference to the clusters.
          const t_Clusters& clusters() const { return clusters_(); }
          //! Reassigns clusters.
          void reassign();
          //! Randomizes both cluster energies and ecis.
          void randomize( typename t_Vector :: value_type _howrandom );

        protected:
          //! Creates the _A and _b matrices for fitting.
          template< class T_MATRIX, class T_VECTOR >
            void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );

          //! The coefficients.
          Indirection::Pointer<t_Matrix> coefficients_;

          //! The classes of equivalent clusters for mixing.
          Indirection::Pointer<t_Clusters> clusters_;

          //! The collapse functor.
          t_Collapse collapse_;

          //! The CE fitting functor.
          Indirection::Pointer<t_CEFit> cefit_;

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
} // namespace LaDa
#include "mixed.impl.h"

#endif
