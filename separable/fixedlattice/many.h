//
//  Version: $Id$
//
#ifndef _CE_MANY_COLLAPSE_H_
#define _CE_MANY_COLLAPSE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>
#include<boost/shared_ptr.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "prepare.h"
#include "colpolicy.h"

namespace LaDa
{
  namespace Traits
  {
    namespace CE
    {
      //! Traits of a collapse functor.
      template< class T_SEPARABLES,
                class T_COLLAPSE,
                class T_MAPPING = ::CE::Mapping::Basic, 
                class T_COEFFICIENTS = boost::numeric::ublas::matrix<types::t_real> >
      struct ManyCollapses 
      {
        //! Type of the configuration matrix.
        typedef T_CONFIGURATIONS t_Configurations;
        //! Type of the Mapping.
        typedef T_MAPPING t_Mapping;
        //! Type of the coefficients.
        typedef typename T_CONFS t_Coefficients;
        //! Type of the individual separable functors.
        typedef typename SeparablesWithMatrixRange< T_SEPARABLES > :: type t_Separables;
        //! Type of the individual collapse functors.
        typedef typename CollapseWithNewSeparables
                         < 
                           T_COLLAPSE, 
                           t_Separables
                         > :: type t_Collapse;
        //! Type of a pair of separables/collapse.
        typedef typename std::pair< t_Collapse, t_Separables > t_FittingPair;
        //! Type of a container of pairs of separables/collapse.
        typedef std::vector< t_FittingPair > t_FittingPairs;

        //! Rebinds type.
        template< class TT_SEPARABLES,
                  class TT_COLLAPSE = T_COLLAPSE , 
                  class TT_MAPPING = T_MAPPING, 
                  class TT_COEFFICIENTS = T_COEFFICIENTS >
          struct rebind
          {
            //! Result type.
            typedef ManyCollapse
                    < 
                      TT_SEPARABLES,
                      TT_COLLAPSE,
                      TT_MAPPING,
                      TT_COEFFICIENTS
                    >  type;
          };
      };

  }

  namespace CE
  {
    //! Collapse functor for fitting CE::Separables  
    template< class T_TRAITS >
      class ManyCollapses
      {
        template< class TT_TRAITS> friend class ManyCollapses;
        public:
          //! Allows rebinding of the collapse function.
          template< class TT_TRAITS > struct rebind
          {
            //! new separable type.
            typedef ManyCollapses< TT_TRAITS > type;
          };
          //! Traits of this functor.
          typedef T_TRAITS t_Traits;
          //! Type of the separable function.
          typedef typename t_Traits :: t_Separables t_Separables;
          //! Type of the mapping function from structures to targets.
          typedef typename t_Traits :: t_Mapping t_Mapping;
          //! Type of the configuration matrix.
          typedef typename t_Traits :: t_Configurations t_Configurations;
          //! Type of the matrices.
          typedef typename t_Separables :: t_Matrix t_Matrix;
          //! Type of the vectors.
          typedef typename t_Separables :: t_Vector t_Vector;

          //! Constructor.
          Collapse() : dim(0), fittingpairs_( new t_FittingPairs ),
                       update_( mapping_ ) {}
          //! Copy Constructor.
          template< class TT_TRAITS>
            ManyCollapses   ( const ManyCollapses< TT_TRAITS > &_c ) 
                     : dim( _c.dim ), fittingpairs_( _c.fittingpairs_ ),
                       mapping_( _c.mapping_ ) {}
          //! Destructor.
          ~ManyCollapses() {}

          //! Creates the fitting matrix and target vector.
          template< class T_MATRIX, class T_VECTOR >
            void operator()( T_MATRIX &_A, T_VECTOR &_b,
                             types::t_unsigned _dim )
              { dim = _dim; create_A_n_b( _A, _b ); regularization()( _A, _b, _dim); }

          //! Evaluates square errors.
          opt::ErrorTuple evaluate() const;
          //! Evaluates square errors for one structure.
          typename t_Matrix :: value_type evaluate( size_t _n ) const;

          //! Updates the scales vector and  normalizes.
          void update_all()
            { foreach( t_FittingPair &pair, *fittingpairs_ ) pair.first.update_all(); }
          //! Updates the scales vector, should update only one dim.
          void update( types::t_unsigned _d );
          //! Does nothing.
          void reset()
            { foreach( t_FittingPair &pair, *fittingpairs_ ) pair.first.reset(); }

          //! Initializes collapse functor.
          template< class T_STRUCTURES >
          void init( const T_STRUCTURES& _strs, const PosToConfs &_postoconfs );
          //! Sets the norm pointer.
          void init( t_Separables& _sep );

          //! Number of configurations.
          size_t nbconfs() const;
          
          //! Returns the number of dimensions.
          size_t dimensions() const;
          //! Returns the number of degrees of liberty (max for all dimension).
          size_t dof() const;
          //! Returns a reference to the separable function;
          t_Separables& separables( size_t _i ) { return (*fittingpairs_)[_i].second; }
          //! Returns a constant reference to the separable function;
          const t_Separables& separables( size_t _i ) const
            { return (*fittingpairs_)[_i].second; }
          //! Returns a reference to the separable function;
          t_Collapses& collapses( size_t _i ) { return (*fittingpairs_)[_i].first; }
          //! Returns a constant reference to the separable function;
          const t_Collapses& collapses( size_t _i ) const
            { return (*fittingpairs_)[_i].first; }
          //! Returns a reference to the mapping.
          t_Mapping& mapping() { return mapping_; }
          //! Returns a reference to the mapping.
          const t_Mapping& mapping() const { return mapping_; }
          //! Returns a reference to the regularization.
          t_RegPolicy& regularization() { return regularization_; }
          //! Returns a constant reference to the regularization.
          const t_RegPolicy& regularization() const { return regularization_; }
          //! Returns a reference to the coefficients.
          typename t_Coefficients::t_Matrix& coefficients()
            { return coefficients(); }
          //! Returns a constant reference to the coefficients.
          const typename t_Coefficients::t_Matrix& coefficients() const
            { return coefficients(); }
          //! Allows manipulation of the coefficients' interface itself.
          typename t_Coefficients&
            coefficients_interface() { return coefficients_; }
          //! Randomizes the coefficients.
          void randomize( typename t_Coefficients :: t_Matrix :: value_type _h)
            { foreach( t_FittingPair &p, *fittingpairs_ ) p.first.randomize( _h ); }
          //! \brief Creates an X vector for fitting, for a single rank.
          //! \details Before the scaling is done, all ranks equal.
          template< class T_VECTOR >
            void create_X( size_t _i, size_t _n, T_VECTOR &_out )
             { dim = _n; create_X( _i, _out ); }
          
          //! Adds a collapse/separable pair and return index.
          size_t addone();

        protected:
          //! Returns the number of degrees of liberty for current dimension.
          size_t current_dof() const;
          //! Creates the _A and _b matrices for fitting.
          template< class T_MATRIX, class T_VECTOR >
            void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );
          //! \brief Creates an X vector for fitting, for a single rank.
          //! \details Before the scaling is done, all ranks equal.
          template< class T_VECTOR >
            void create_X( size_t _i, T_VECTOR &_out );


          //! Holds current dimension being fitted.
          size_t dim;
          //! Pointer to separable function being minimized.
          boost::shared_ptr< t_Fittingpairs > fittingpairs_;
          //! \brief The mapping from target values to symetrically equivalent
          //!        structures.
          t_Mapping mapping_;
          //! Update policy.
          t_UpdatePolicy update_;
          //! Regularization.
          t_RegPolicy regularization_;
      };
   
    //! Prints description of the separables function in the collapse functor.
    template<class T_TRAITS>
      std::ostream& operator<<( std::ostream& _stream,
                                const ManyCollapses<T_TRAITS>& _col );

    //! Saves the state of a collapse object.
    class CollapseState
    {
      public:
        template< class T_COLLAPSE >
          void operator=( const T_COLLAPSE& _c )
          {
            coefficients_ = _c.coefficients();
            norms_ = _c.separables().norms;
          }
        template< class T_COLLAPSE >
          void reset( T_COLLAPSE& _c ) const
          {
            _c.coefficients() = coefficients_;
            _c.separables().norms = norms_;
          }

      protected:
        //! Coefficients to save.
        typedef boost::numeric::ublas::matrix<types::t_real> t_Coefficients;
        //! Norms to save.
        typedef boost::numeric::ublas::vector<types::t_real> t_Norms;
        //! Coefficients to save.
        t_Coefficients coefficients_;
        //! Norms to save.
        t_Norms norms_;
    };

  }
} // namespace LaDa
#include "collapse.impl.h"

#endif
