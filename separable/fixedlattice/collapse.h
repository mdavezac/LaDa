//
//  Version: $Id$
//
#ifndef _CE_COLLAPSE_H_
#define _CE_COLLAPSE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "prepare.h"

namespace CE
{
  //! \cond
  namespace Methods
  {
    template< class T_COLLAPSE, class T_SEPARABLES, class T_STRUCTURES >
      opt::ErrorTuple check_one( const T_SEPARABLES &_separables,
                                 const T_COLLAPSE &_collapse,
                                 const T_STRUCTURES &_str,
                                 size_t _n, bool _verbose = false );
  }
  namespace Policy
  {
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      class LowMemUpdate;
  }
  //! \endcond
} // end of CE namespace

namespace Traits
{
  namespace CE
  {
    //! Traits of a collapse functor.
    template< class T_SEPARABLES,
              class T_MAPPING = ::CE::Mapping::SymEquiv, 
              template<class, class, class>
                class T_UPDATEPOLICY = ::CE::Policy::LowMemUpdate >
    struct Collapse 
    {
      //! Type of the configuration matrix.
      typedef boost::numeric::ublas::matrix<size_t> t_iMatrix;
      //! Type of the Mapping.
      typedef T_SEPARABLES t_Separables;
      //! Type of the Mapping.
      typedef T_MAPPING t_Mapping;
      //! Type of the Policy.
      typedef T_UPDATEPOLICY<t_Separables, t_Mapping, t_iMatrix> t_UpdatePolicy;
    };
  }
} // end of traits namespace.

namespace CE
{
  //! Collapse functor for fitting CE::Separables  
  template< class T_TRAITS >
    class Collapse
    {
      template< class T_COLLAPSE, class TT_SEPARABLES, class T_STRUCTURES > friend
        opt::ErrorTuple Methods::check_one( const TT_SEPARABLES &_separables,
                                            const T_COLLAPSE &_collapse,
                                            const T_STRUCTURES &_str,
                                            size_t _n, bool _verbose );
      public:
        //! Traits of this functor.
        typedef T_TRAITS t_Traits;
        //! Type of the separable function.
        typedef typename t_Traits :: t_Separables t_Separables;
        //! Type of the mapping function from structures to targets.
        typedef typename t_Traits :: t_Mapping t_Mapping;
        //! Type of the configuration matrix.
        typedef typename t_Traits :: t_iMatrix t_iMatrix;
        //! Type of the update policy.
        typedef typename t_Traits :: t_UpdatePolicy t_UpdatePolicy;
        //! Type of the matrices.
        typedef typename t_Separables :: t_Matrix t_Matrix;
        //! Type of the vectors.
        typedef typename t_Separables :: t_Vector t_Vector;

        //! \brief The mapping from target values to symetrically equivalent
        //!        structures.
        t_Mapping mapping;

        //! Constructor.
        Collapse() : dim(0), separables_(NULL),
                     update_( mapping, configurations_ ) {}
        //! Destructor.
        ~Collapse() {}

        //! Creates the fitting matrix and target vector.
        void operator()( t_Matrix &_A, t_Vector &_b,
                         types::t_unsigned _dim )
          { dim = _dim; create_A_n_b( _A, _b ); }
        //! Evaluates square errors.
        opt::ErrorTuple evaluate();

        //! Updates the scales vector and  normalizes.
        void update_all();
        //! Updates the scales vector, should update only one dim.
        void update( types::t_unsigned _d );
        //! Does nothing.
        void reset() {}

        //! Initializes collapse functor.
        template< class T_STRUCTURES >
        void init( const T_STRUCTURES& _strs, const PosToConfs &_postoconfs );
        //! Sets the norm pointer.
        void init( t_Separables& _sep );

        //! Reference to configuration matrix.
        const t_iMatrix& configurations() const { return configurations_; }
        
        //! Returns the number of dimensions.
        size_t dimensions() const { return separables_->dimensions(); }
        //! Returns the number of degrees of liberty (per dimension).
        size_t dof() const { return separables_->dof(); }
        //! Returns a constant reference to the separable function;
        t_Separables& separables() { return *separables_; }

      protected:
        //! Creates the _A and _b matrices for fitting.
        void create_A_n_b( t_Matrix &_A, t_Vector &_b );
        //! Finds scaling factor for that conf, collapsed dimension, and rank.
        typename t_Vector::value_type factor( size_t _kv, size_t _r, size_t _d );
        //! \brief Creates an X vector for fitting, for a single rank.
        //! \details Before the scaling is done, all ranks equal.
        void create_X( size_t _i, t_Vector &_out );


        //! The configurations, arranged in columns.
        t_iMatrix configurations_;
        //! Holds current dimension being fitted.
        size_t dim;
        //! holds the sep. func. split up by rank and confs.
        t_Matrix scales;
        //! Pointer to separable function being minimized.
        t_Separables *separables_;
        //! Update policy.
        t_UpdatePolicy update_;
    };

  namespace Policy
  {
    //! \brief Low Memory update scheme.
    //! \details Only the values of each rank of the separable function is kept
    //!          for each configuration. As such, updates happen over all
    //!          dimensions always.
    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      class LowMemUpdate
      {
        public:
          //! Type of the separable functions.
          typedef T_SEPARABLES t_Separables;
          //! Type of the mapping from configurations to configurations.
          typedef T_MAPPING t_Mapping;
          //! Type of the mapping from configurations to configurations.
          typedef T_CONFS t_iMatrix;
          //! Type of the vectors.
          typedef typename t_Separables :: t_Vector t_Vector;
          //! Type of the matrices.
          typedef typename t_Separables :: t_Matrix t_Matrix;

          //! the constructor.
          LowMemUpdate   ( const t_Mapping &_mapping, const t_iMatrix &_confs )
                       : mapping_(_mapping), configurations_(_confs) {}
          //! Destructor.
          ~LowMemUpdate() {}

          //! Updates all dimension.
          void operator()();
          //! Updates only dimension \a _dim. 
          void operator()( size_t _dim ) { operator()(); }
          //! Returns the factor for configurations _kv, and rank _r, and dimension.
          typename t_Vector :: value_type factor( size_t _kv, size_t _r, size_t _dim) const;
          //! Updates pointer to separable function.
          void init( const t_Separables& _sep );

        protected:
          //! Recomputes factor from nothing.
          typename t_Vector :: value_type
            factor_from_scratch( size_t _kv, size_t _r, size_t _dim) const;

          //! Holds the sep. func. split up by rank and confs.
          t_Matrix scales_;
          //! A reference to the config. mapping.
          const t_Mapping &mapping_;
          //! A reference to the config. mapping.
          const t_iMatrix &configurations_;
          //! A reference to the separable function.
          const t_Separables * separables_;
      };

    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      class HighMemUpdate : protected LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS >
      {
        public:
          //! Type if the base class.
          typedef LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS > t_Base;
          //! Type of the separable functions.
          typedef T_SEPARABLES t_Separables;
          //! Type of the mapping from configurations to configurations.
          typedef T_MAPPING t_Mapping;
          //! Type of the mapping from configurations to configurations.
          typedef T_CONFS t_iMatrix;
          //! Type of the vectors.
          typedef typename t_Separables :: t_Vector t_Vector;
          //! Type of the matrices.
          typedef typename t_Separables :: t_Matrix t_Matrix;

          //! Constructor.
          HighMemUpdate   ( const t_Mapping &_mapping, const t_iMatrix &_confs )
                        : t_Base( _mapping, _confs ) {}
          //! Updates all dimension.
          void operator()();
          //! Updates only dimension \a _dim. 
          void operator()( size_t _dim );
          //! Returns the factor for configurations _kv, and rank _r, and dimension.
          typename t_Vector :: value_type factor( size_t _kv, size_t _r, size_t _dim) const;
          //! Updates pointer to separable function.
          void init( const t_Separables& _sep );

        protected:
          //! Updates scales_ from dimsplit.
          void update_scales();
          //! Returns the factor for configurations _kv, and rank _r, and dimension.
          typename t_Vector :: value_type
            factor_from_scratch( size_t _kv, size_t _r, size_t _dim) const;

          //! Holds the sep. func. split along confs, ranks, and dimensions.
          std::vector< t_Matrix > dimsplit_;
          using t_Base :: scales_;
          using t_Base :: mapping_;
          using t_Base :: configurations_;
          using t_Base :: separables_;
      };
  } // end of Policy namespace.

}

#include "collapse.impl.h"

#endif
