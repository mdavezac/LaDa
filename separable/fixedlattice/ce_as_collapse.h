//
//  Version: $Id$
//
#ifndef _CE_MANY_H_
#define _CE_MANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "colpolicy.h"

namespace CE
{
  //! \cond
  namespace details
  {
    template < class T_FIT > Regularization; 
  }
  //! \endcond

  //! Wraps a CE::Fit with the interface of a collapse functor.
  template< class T_POLICY >
  class CEasCollapse 
  {
    public:
      //! Type of the fitting functor.
      typedef Fit< T_POLICY > t_Fit;
      //! Type of the mapping.
      typedef Mapping::Basic t_Mapping;
      //! Type of the mapping.
      typedef details::Regularization< t_Fit > t_Regularization;
     
      //! constructor.
      CEasCollapse( CEasCollapse& _c ) : fit_( _c.fit_ ) 
        { clusters_.reset( new t_Clusters ); }
      //! constructor.
      CEasCollapse( t_Fit& _fit ) : fit_( _fit ), clusters_( _c.clusters_ ) {}

      //! Evaluates square errors for one structure.
      typename t_Matrix :: value_type evaluate( size_t _n ) const;

      //! Updates the scales vector and  normalizes.
      void update_all() {}
      //! Updates the scales vector, should update only one dim.
      void update( types::t_unsigned _d ) {}
      //! Does nothing.
      void reset() {}

      //! Number of configurations.
      size_t nbconfs() const { return 1; }
      //! Returns the number of dimensions.
      size_t dimensions() const { return 1; }
      //! Returns the number of degrees of liberty (per dimension).
      size_t dof() const { return clusters_->size(); }

      //! Returns a reference to the mapping.
      t_Mapping& mapping() { return mapping_; }
      //! Returns a reference to the mapping.
      const t_Mapping& mapping() const { return mapping_; }
      //! Returns a reference to the regularization.
      t_Regularization& regularization() { return regularization_; }
      //! Returns a constant reference to the regularization.
      const t_Regularization& regularization() const { return regularization_; }
      //! Randomizes the coefficients.
      void randomize( types::t_real _howrandom );
      //! Creates an X vector for fitting, for a single rank.
      template< class T_VECTOR > void create_X( size_t _i, size_t, T_VECTOR &_out )
        { std::copy( cefit.pis[_i].begin(), cefit.pis[_i].end(), X.begin() ); }

    protected:
      //! Reference to the CE::Fit object.
      t_Fit &fit_;
      //! The classes of equivalent clusters for mixing.
      boost::shared_ptr<t_Clusters> clusters_;
      //! the mapping from configuration to clusters.
      t_Mapping mapping_;
      //! The regularization.
      t_Regularization regularization_;
  };


  template< class T_POLICY >
    void CEasCollapse<T_POLICY> :: randomize( types::t_real _howrandom )
    {
    }

}
#endif

