//
//  Version: $Id$
//
#ifndef _CE_MANY_H_
#define _CE_MANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "colpolicy.h"
#include <opt/random.h>

namespace Traits
{
  namespace CE
  {
    //! A traits class for wrapping a CE fitting functor as a collapse functor.
    template< class T_POLICY,
              template<class> class T_CEFIT,
              class T_COEFFICIENTS >
    struct CEasCollapse
    {
      //! Type of the policy.
      typedef T_POLICY t_Policy;
      //! Type of the fitting functor.
      typedef T_CEFIT<t_Policy> t_CEFit;
      //! Type of the coefficients.
      typedef T_COEFFICIENTS t_Coefficients;
      //! Helps to obtain a differently typed fitting functor.
      template< class TT_POLICY > struct rebindCEFit
      { 
        //! Type of the regularization policy.
        typedef T_CEFIT< TT_POLICY > other;
      };
    };
  }
}

namespace CE
{
  //! \cond
  namespace details
  {
    template< class T_CEFIT > class Regularization;
  }
  //! \endcond

  //! Wraps a CE::Fit with the interface of a collapse functor.
  template< class T_TRAITS >
  class CEasCollapse 
  {
    public:
      //! Type of the traits
      typedef T_TRAITS t_Traits;
      //! Type of the policy.
      typedef typename t_Traits :: t_Policy t_Policy;
      //! Type of the fitting functor.
      typedef typename t_Traits :: t_CEFit t_Fit;
      //! Type of the mapping.
      typedef Mapping::Basic t_Mapping;
      //! Type of the mapping.
      typedef details::Regularization< t_CEFit > t_Regularization;
      //! Type of the coefficients.
      typedef t_Traits :: t_Coefficients t_Coefficients;
      //! Type of the matrices.
      typedef t_Coefficients :: t_Matrix t_Matrix;
     
      //! constructor.
      CEasCollapse() { clusters_.reset( new t_Clusters ); }
      //! Copy constructor.
      CEasCollapse   ( CEasCollapse& _c )
                   : fit_( _c.fit_ ), clusters_(_c.clusters_),
                     mapping_( _c.mapping_), regularization( _c.regularization_ ),
                     coefficients_(_c.coefficients_) {}
      //! constructor.
      CEasCollapse( t_Fit& _fit ) : fit_( _fit ), clusters_( _c.clusters_ ) {}

      //! Evaluates square errors for one structure.
      typename t_Matrix::value_type evaluate( size_t _n ) const
       { return bblas::inner_prod( coefficients(), cefit().pis[ _n ] ); }

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
      //! Creates an X vector for fitting, for a single rank.
      template< class T_VECTOR > void create_X( size_t _i, size_t, T_VECTOR &_out )
        { std::copy( cefit.pis[_i].begin(), cefit.pis[_i].end(), X.begin() ); }
      //! Reassigns clusters.
      void reassign();
      //! Returns a reference to the coefficients.
      t_Matrix& coefficients() { return coefficients_(); }
      //! Returns a constant reference to the coefficients.
      const t_Matrix& coefficients() const { return coefficients_(); }
      //! Allows manipulation of the coefficients' interface itself.
      t_Coefficients& coefficients_interface() { return coefficients_; }
      //! Randomizes the coefficients.
      void randomize( typename t_Matrix :: value_type _howrandom )

    protected:
      //! Reference to the CE::Fit object.
      t_Fit &fit_;
      //! The classes of equivalent clusters for mixing.
      boost::shared_ptr<t_Clusters> clusters_;
      //! the mapping from configuration to clusters.
      t_Mapping mapping_;
      //! The regularization.
      t_Regularization regularization_;
      //! The coefficients.
      t_Coefficients coefficients_;

  };

  template< class T_TRAITS > void CEasCollapse<T_TRAITS> :: reassign()
  {
    namespace bl = boost::lambda;
    __DEBUGTRYBEGIN
    __ASSERT( coefficients().size() != clusters_->size(), "Inconsistent sizes.\n")
 
    typename T_VECTOR :: const_iterator i_eci = coefficients().begin();
    foreach( t_Clusters::value_type & _clusters, clusters() )
    {
      foreach( ::CE::Cluster & _cluster, _clusters ) _cluster.eci = *i_eci;
      ++i_eci;
    }
    __DEBUGTRYEND(,"Error in CEasCollapse::reassign().\n" )
  }
  template< class T_TRAITS >
    void CEasCollapse<T_TRAITS> :: randomize( typename t_Matrix :: value_type _howrandom )
    {
      typename t_Matrix :: iterator i_c = coefficients().begin();
      typename t_Matrix :: iterator i_c_end = coefficients().end();
      for(; i_c != i_c_end; ++i_c )
        *i_c = t_Type( opt::random::rng() - 5e-1 ) * _howrandom;
    }

  namespace details
  {
    //! \brief Wrapper around CE::RegulatedFit regularization.
    template< class T_CEFIT >
    class Regularization : public CE::Policy::NoReg< T_CEFIT >
    {
      public:
        //! Helps to obtain a differently typed regularization.
        template< class TT_CEFIT > struct rebind
        { 
          //! Type of the regularization policy.
          typedef Regularization< TT_CEFIT > other;
        };
        //! Type of the separable function.
        typedef T_CEFIT t_CEFit;

        //! Constructor
        Regularization() : CE::Policy::NoReg<t_CEFit>() {};
        //! Copy Constructor.
        template <class TT_CEFIT> 
          Regularization   ( const Regularization<TT_CEFIT> &_c )
                         : CE::Policy::NoReg<t_CEFit>( _c ) {};
        //! Copy Constructor.
        template <class TT_SEPARABLES> 
          Regularization   ( const CE::Policy::NoReg<TT_CEFIT> &_c )
                         : CE::Policy::NoReg<t_CEFit>( _c ) {};

        //! Would modify A matrix and b vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &_m, T_VECTOR &_v, size_t _dim )
            { if( _dim == 1 ) separables_.other_A_n_b( _m, _v ); };

        using CE::Policy::NoReg<t_CEFit> :: init;
      protected:
        using CE::Policy::NoReg<t_CEFit> :: separables_;
    };
  }


}
#endif

