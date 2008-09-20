//
//  Version: $Id$
//
#ifndef _SEPARABLES_FIXEDLATTICE_CE_AS_COLLAPSE_H_
#define _SEPARABLES_FIXEDLATTICE_CE_AS_COLLAPSE_H_

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
    template< class T_CEFIT, class T_COEFFICIENTS >
    struct CEasCollapse
    {
      //! Type of the fitting functor.
      typedef T_CEFIT t_CEFit;
      //! Type of the policy.
      typedef typename t_CEFit :: t_Policy t_Policy;
      //! Type of the coefficients.
      typedef T_COEFFICIENTS t_Coefficients;
      //! Helps to obtain a differently typed class.
      template< class TT_CEFIT, class TT_COEFFICIENTS > struct rebind
      { 
        //! Resulting type.
        typedef CEasCollapse< TT_CEFIT, TT_COEFFICIENTS > type;
      };
      //! Fake nullary functor
      template< class T > struct rebind_with_new_separables
      {
        //! The same type.
        typedef CEasCollapse< T_CEFIT, T_COEFFICIENTS > type;
      };
    };
    //! \cond
    namespace details
    {
      template< class T_MAPPING, class T_POLICY, 
                class T_COEFFICIENTS, class T_VECTOR,
                class T_FIT > struct CEasSeparables;
    }
    //! \endcond
    //! Metafunction for wrapping a CE traits as Separables function traits.
    template< class T_FIT,
              class T_COEFFICIENTS = ::CE::Policy::MatrixRangeCoefficients >
      struct WrapCEasSeparables
      {
        protected:
          //! Type of the fitting function.
          typedef T_FIT t_Fit;
          //! Type of the coefficients.
          typedef T_COEFFICIENTS t_Coefficients;
          //! Type of the policy
          typedef typename t_Fit :: t_Policy t_Policy;
          //! A fake mapping with the correct types and constants.
          struct t_Mapping 
          {
            //! The degrees of liberty per "rank".
            const static size_t D = 1;
          };
        public:
          //! Result type.
          typedef details::CEasSeparables
                  < 
                    t_Mapping, t_Policy, t_Coefficients, 
                    boost::numeric::ublas::vector
                    < 
                      typename t_Coefficients::t_Matrix::value_type
                    >,
                   T_FIT 
                  > type;
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
      //! Helps to obtain a differently typed class.
      template< class TT_TRAITS > struct rebind
      { 
        //! Resulting type.
        typedef CEasCollapse< TT_TRAITS > type;
      };
      //! Type of the traits
      typedef T_TRAITS t_Traits;
      //! Type of the policy.
      typedef typename t_Traits :: t_Policy t_Policy;
      //! Type of the fitting functor.
      typedef typename t_Traits :: t_CEFit t_CEFit;
      //! Type of the mapping.
      typedef Mapping::Basic t_Mapping;
      //! Type of the mapping.
      typedef details::Regularization< t_CEFit > t_Regularization;
      //! Type of the coefficients.
      typedef typename t_Traits :: t_Coefficients t_Coefficients;
      //! Type of the matrices.
      typedef typename t_Coefficients :: t_Matrix t_Matrix;
      //! The type of the clusters.
      typedef typename t_CEFit :: t_Clusters t_Clusters;
     
      //! constructor.
      CEasCollapse() : clusters_( new t_Clusters ), cefit_(NULL) {}
      //! constructor.
      CEasCollapse   ( t_CEFit& _fit ) 
                   : cefit_( _fit ), clusters_( new t_Clusters ) {}
      //! Copy constructor.
      CEasCollapse   ( const CEasCollapse& _c )
                   : cefit_( _c.cefit_ ), clusters_(_c.clusters_),
                     mapping_( _c.mapping_), regularization_( _c.regularization_ ) {}
      
      //! Sets the pointer to the fitting class.
      void init( t_CEFit &_fit ) { cefit_ = &_fit; }

      //! Evaluates square errors for one structure.
      typename t_Matrix::value_type evaluate( size_t _n ) const;

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
        { std::copy( cefit.pis[_i].begin(), cefit.pis[_i].end(), _out.begin() ); }
      //! Reassigns clusters.
      void reassign();
      //! Returns a reference to the coefficients.
      t_CEFit& separables() { return cefit(); }
      //! Returns a reference to the coefficients.
      const t_CEFit& separables() const { return cefit(); }
      //! Returns a reference to the coefficients.
      t_Matrix& coefficients() { return separables().coefficients(); }
      //! Returns a constant reference to the coefficients.
      const t_Matrix& coefficients() const { return separables().coefficients(); }
      //! Allows manipulation of the coefficients' interface itself.
      t_Coefficients& coefficients_interface() { return separables().coefficients_interface(); }
      //! Randomizes the coefficients.
      void randomize( typename t_Matrix :: value_type _howrandom );
      //! Initializes clusters and mapping.
      template< class T_CLUSTERS > void init( const T_CLUSTERS &_Clusters );
      //! Returns a constant reference to the clusters.
      const t_Clusters& clusters() const { return *clusters_; }

    protected:
      //! Returns a reference to the fitting object.
      t_CEFit& cefit()
      { 
        __ASSERT( cefit_ == NULL, "Null pointer.\n" )
        return *cefit_; 
      }
      //! Returns a constant reference to the fitting object.
      const t_CEFit& cefit() const
      { 
        __ASSERT( cefit_ == NULL, "Null pointer.\n" )
        return *cefit_; 
      }
      //! Reference to the CE::Fit object.
      t_CEFit *cefit_;
      //! The classes of equivalent clusters for mixing.
      boost::shared_ptr<t_Clusters> clusters_;
      //! the mapping from configuration to clusters.
      t_Mapping mapping_;
      //! The regularization.
      t_Regularization regularization_;

  };

  //! Prints description of the clusters in the collapse functor.
  template<class T_TRAITS>
  std::ostream& operator<<( std::ostream& _stream, const CEasCollapse<T_TRAITS>& _col )
  {
    typedef typename T_TRAITS :: t_CEFit :: t_Clusters t_Clusters;
    foreach( const typename t_Clusters :: value_type &cluster_class, _col.clusters() )
      std::cout << cluster_class.front() << " D=" << cluster_class.size() << "\n";
    return _stream;
  }

  //! \brief Wraps a CE::Fit class to include separables style metafunctions.
  template< class T_CEFIT, 
            class T_TRAITS = typename Traits::CE::WrapCEasSeparables<T_CEFIT>::type >
    class CEasSeparables : public T_CEFIT
    {
      public:
        //! Type of the traits.
        typedef T_TRAITS t_Traits;
        //! Type of the fitting function.
        typedef T_CEFIT t_CEFit;
        //! Allows rebinding of the fake separables function traits.
        template< class TT_TRAITS >
          struct rebind
          {
            //! new type.
            typedef CEasSeparables< T_CEFIT, TT_TRAITS > type;
          };
        //! Type of the coefficients.
        typedef typename t_Traits :: t_Coefficients t_Coefficients;

        
        //! Constructor.
        CEasSeparables() {}
        //! Copy Constructor.
        CEasSeparables   ( const CEasSeparables & _c ) 
                       : T_CEFIT( _c ), coefficients_(_c.coefficients_) {}
        //! Returns a reference to the coefficients.
        typename t_Coefficients :: t_Matrix & coefficients()
          { return coefficients_(); }
        //! Returns a constant reference to the coefficients.
        const typename t_Coefficients :: t_Matrix & coefficients() const
          { return coefficients_(); }
        //! Allows manipulation of the coefficients' interface itself.
        t_Coefficients& coefficients_interface() { return coefficients_; }

      protected:
         //! The coefficients.
         t_Coefficients coefficients_;
    };

}

#include "ce_as_collapse.impl.h"
#endif

