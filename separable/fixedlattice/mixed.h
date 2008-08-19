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
    template< class T_COLTRAITS, class T_CEBASE >
     struct MixedApproach
     {
       //! Original separable traits.
       typedef typename T_COLTRAITS :: t_Separables :: t_Traits t_OrigSepTraits;
       //! Separable function traits.
       typedef Separables
               < 
                 typename t_OrigSepTraits :: t_Mapping,
                 typename t_OrigSepTraits :: t_Policy,
                 ::CE::Policy::MatrixRangeCoefficients,
                 typename t_OrigSepTraits :: t_Vector
               > 
               t_SepTraits;
       //! Type of the separable function.
       typedef ::CE::Separables< t_SepTraits > t_Separables;
       //! Type of the configuration matrix.
       typedef typename T_COLTRAITS :: t_Configurations t_Configurations;
       //! Type of the Mapping.
       typedef typename T_COLTRAITS :: t_Mapping t_Mapping;
       //! Type of the Regulations Policy
       typedef typename T_COLTRAITS :: t_RegPolicy
                                    ::template rebind< t_Separables > :: other 
          t_RegPolicy;
       //! Type of the Policy.
       typedef typename T_COLTRAITS :: t_UpdatePolicy
                           ::template rebind< t_Separables,
                                           typename T_COLTRAITS :: t_Mapping, 
                                           typename T_COLTRAITS :: t_Configurations >
                           :: other t_UpdatePolicy;
       //! collapse traits.
       typedef Collapse
               < 
                 t_Separables,
                 t_Mapping,
                 t_RegPolicy,
                 t_Configurations,
                 t_UpdatePolicy
               > 
               t_ColTraits;
       //! Type of the collapse functor base.
       typedef ::CE::Collapse< t_ColTraits > t_Collapse;
       //! CE fit base.
       typedef T_CEBASE t_CEBase;
     };

  } // end of CE namespace 

} // end of Traits namespace
namespace CE
{

  template< class T_TRAITS >
    class MixedApproach : protected T_TRAITS :: t_Collapse,
                          protected T_TRAITS :: t_CEBase
    {
      template< class TT_TRAITS > friend class MixedApproach;
      public:
        //! Type of the traits.
        typedef T_TRAITS t_Traits;
        //! Type of the separable collapse functor.
        typedef typename t_Traits :: t_Collapse t_ColBase;
        //! Type of the CE fitting base class.
        typedef typename t_Traits :: t_CEBase  t_CEBase;
        //! Type of the separable function.
        typedef typename t_Traits :: t_Separables  t_Separables;
        //! Type of the configuration matricex.
        typedef typename t_Traits :: t_ColTraits :: t_Configurations t_Configurations;
        //! \brief Type of the matrix range.
        //! \details Necessary interface for minimizer.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Matrix t_Matrix;
        //! Type of the vectors.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Vector t_Vector;
        //! Type of the classes of equivalent clusters.
        typedef std::vector< std::vector< ::CE::Cluster > > t_Clusters;


        //! Constructor.
        MixedApproach() : t_ColBase(), t_CEBase() 
          { clusters_.reset( new t_Clusters ); coefficients_.reset( new t_Matrix ); }
        //! Constructor.
        template< class TT_TRAITS >
          MixedApproach   ( const MixedApproach< TT_TRAITS > &_c )
                        : t_ColBase( _c ), t_CEBase( _c ),
                          separables_( _c.separables_ ), 
                          clusters_( _c.clusters_ ), 
                          coefficients_( _c.coefficients_ ) 
            { init( _c.separables().ranks(), _c.separables().dimensions() ); }
        //! Destructor.
        ~MixedApproach() {}

        //! Returns a reference to the CE fitting base.
        t_CEBase& CEFit() { return *( (t_CEBase*) this ); }
        //! Returns a constant reference to the CE fitting base.
        const t_CEBase& CEFit() const { return *( (t_CEBase*) this ); }
        //! Returns a reference to the CE fitting base.
        t_ColBase& Collapse() { return *( (t_ColBase*) this ); }
        //! Returns a constant reference to the CE fitting base.
        const t_ColBase& Collapse() const { return *( (t_ColBase*) this ); }


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

        //! Reference to configuration matrix.
        const t_Configurations& configurations() const
          { return t_ColBase::configurations(); }
        
        //! Returns the number of dimensions.
        size_t dimensions() const { return t_ColBase::dimensions(); }
        //! Returns the number of degrees of liberty (per dimension).
        size_t dof() const { return t_ColBase::dof() + t_CEBase::dof(); }
       
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
        typename t_ColBase::t_Mapping& mapping() { return t_ColBase::mapping(); }
        //! Returns a reference to the mapping.
        const typename t_ColBase::t_Mapping& mapping() const
          { return t_ColBase::mapping(); }
        //! Returns a reference to the regularization.
        typename t_ColBase::t_RegPolicy& regularization() 
          { return t_ColBase::regularization(); }
        //! Returns a constant reference to the regularization.
        const typename t_ColBase::t_RegPolicy& regularization() const
          { return t_ColBase::regularization(); }
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
    };

  //! Prints mixed-approach description to a stream.
  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream, const MixedApproach<T_TRAITS> &_col );

  //! \brief Removes cluster classes for wich at least one custer is contained
  //!        within \a _positions and the origin.
  template< class T_POSITIONS, class T_CLUSTERS >
    void remove_contained_clusters( const T_POSITIONS &_positions, T_CLUSTERS &_clusters );
} // end of CE namspace

#include "mixed.impl.h"

#endif
