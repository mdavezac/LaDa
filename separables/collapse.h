#ifndef _CE_COLLAPSE_H_
#define _CE_COLLAPSE_H_

#include "LaDaConfig.h"

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include <misc/types.h>
#include <opt/debug.h>
#include <opt/indirection.h>

#include "prepare.h"
#include "colpolicy.h"

namespace LaDa
{
  //! \cond
  namespace CE
  {
    namespace Methods
    {
      template< class T_COLLAPSE, class T_SEPARABLES, class T_STRUCTURES >
        opt::ErrorTuple check_one( const T_SEPARABLES &_separables,
                                   const T_COLLAPSE &_collapse,
                                   const T_STRUCTURES &_str,
                                   size_t _n, bool _verbose = false );
    }
  } // end of CE namespace
  //! \endcond

  namespace Traits
  {
    namespace CE
    {
      //! Traits of a collapse functor.
      template< class T_SEPARABLES,
                class T_MAPPING = ::LaDa::CE::Mapping::SymEquiv, 
                class T_REGULARIZATIONPOLICY = ::LaDa::CE::Policy::NoReg<T_SEPARABLES>,
                class T_CONFS = boost::numeric::ublas::matrix<size_t>,
                class T_UPDATEPOLICY
                  = ::LaDa::CE::Policy::LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS > >
      struct Collapse 
      {
        //! Type of the configuration matrix.
        typedef T_CONFS t_Configurations;
        //! Type of the Mapping.
        typedef T_SEPARABLES t_Separables;
        //! Type of the coefficients.
        typedef typename t_Separables :: t_Coefficients :: t_Matrix t_Coefficients;
        //! Type of the Mapping.
        typedef T_MAPPING t_Mapping;
        //! Type of the Regulations Policy
        typedef T_REGULARIZATIONPOLICY t_RegPolicy;
        //! Type of the Policy.
        typedef T_UPDATEPOLICY t_UpdatePolicy;

        //! Rebinds type.
        template< class TT_SEPARABLES = t_Separables,
                  class TT_MAPPING = t_Mapping,
                  class TT_REGULARIZATIONPOLICY = t_RegPolicy,
                  class TT_CONFS = t_Configurations, 
                  class TT_UPDATEPOLICY = t_UpdatePolicy > 
          struct rebind
          {
            //! Result type.
            typedef Collapse< TT_SEPARABLES, TT_MAPPING,
                              TT_REGULARIZATIONPOLICY, TT_CONFS, 
                              TT_UPDATEPOLICY >  type;
          };
        //! Rebinds the traits with new separables.
        template< class T_NEWSEP > struct rebind_with_new_separables;
      };

      //! Rebinds collapse with new separables.
      template< class T_COLLAPSE, class T_NEWSEP > struct CollapseWithNewSeparables
      { 
        protected:
          //! The new traits.
          typedef typename T_COLLAPSE::t_Traits::template rebind_with_new_separables
                           < 
                             T_NEWSEP 
                           > :: type t_NewTraits;
        public:
          //! Result type.
          typedef typename T_COLLAPSE :: template rebind
                           < 
                             t_NewTraits
                           > :: type type;
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
        template< class TT_TRAITS> friend class Collapse;
        public:
          //! Allows rebinding of the collapse function.
          template< class TT_TRAITS > struct rebind
          {
            //! new separable type.
            typedef Collapse< TT_TRAITS > type;
          };
          //! Traits of this functor.
          typedef T_TRAITS t_Traits;
          //! Type of the separable function.
          typedef typename t_Traits :: t_Separables t_Separables;
          //! Type of the mapping function from structures to targets.
          typedef typename t_Traits :: t_Mapping t_Mapping;
          //! Type of the configuration matrix.
          typedef typename t_Traits :: t_Configurations t_Configurations;
          //! Type of the update policy.
          typedef typename t_Traits :: t_UpdatePolicy t_UpdatePolicy;
          //! Type of the update policy.
          typedef typename t_Traits :: t_RegPolicy t_RegPolicy;
          //! Type of the matrices.
          typedef typename t_Separables :: t_Matrix t_Matrix;
          //! Type of the vectors.
          typedef typename t_Separables :: t_Vector t_Vector;

          //! Constructor.
          Collapse() : dim(0), configurations_( new t_Configurations ) {}
          //! Copy Constructor.
          Collapse   ( const Collapse &_c ) 
                   : dim( _c.dim ), separables_( _c.separables_ ),
                     configurations_( _c.configurations_ ),
                     mapping_( _c.mapping_ ) {}
          //! Destructor.
          ~Collapse() {}
          //! Assignement operator.
          template< class TT_TRAITS >
            void operator=( const Collapse<TT_TRAITS> &_c );

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
          void update_all();
          //! Updates the scales vector, should update only one dim.
          void update( types::t_unsigned _d );
          //! Does nothing.
          void reset() {}

          //! Initializes collapse functor.
          template< class T_STRUCTURES >
          void init( const T_STRUCTURES& _strs, const PosToConfs &_postoconfs );
          //! Sets the norm pointer.
  //       void init( t_Separables& _sep );

          //! Number of configurations.
          size_t nbconfs() const { return configurations().size2(); }
          
          //! Returns the number of dimensions.
          size_t dimensions() const { return separables().dimensions(); }
          //! Returns the number of degrees of liberty (per dimension).
          size_t dof() const { return separables().dof(); }
          //! Returns a reference to the separable function;
          t_Separables& separables() { return separables_(); }
          //! Returns a constant reference to the separable function;
          const t_Separables& separables() const { return separables_(); }
          //! Returns a reference to the mapping.
          t_Mapping& mapping() { return mapping_; }
          //! Returns a reference to the mapping.
          const t_Mapping& mapping() const { return mapping_; }
          //! Returns a reference to the regularization.
          t_RegPolicy& regularization() { return regularization_; }
          //! Returns a constant reference to the regularization.
          const t_RegPolicy& regularization() const { return regularization_; }
          //! Returns a reference to the coefficients.
          typename t_Separables::t_Coefficients::t_Matrix& coefficients()
            { return separables().coefficients(); }
          //! Returns a constant reference to the coefficients.
          const typename t_Separables::t_Coefficients::t_Matrix& coefficients() const
            { return separables().coefficients(); }
          //! Allows manipulation of the coefficients' interface itself.
          typename t_Separables :: t_Coefficients&
            coefficients_interface() { return separables().coefficients_interface(); }
          //! Randomizes the coefficients.
          void randomize( typename t_Vector :: value_type _howrandom )
            { separables().randomize( _howrandom ); }
          //! \brief Creates an X vector for fitting, for a single rank.
          //! \details Before the scaling is done, all ranks equal.
          template< class T_VECTOR >
            void create_X( size_t _i, size_t _n, T_VECTOR &_out )
             { dim = _n; create_X( _i, _out ); }

        protected:
          //! Reference to configuration matrix.
          t_Configurations& configurations() { return *configurations_; }
          //! Constant reference to configuration matrix.
          const t_Configurations& configurations() const { return *configurations_; }
          //! Creates the _A and _b matrices for fitting.
          template< class T_MATRIX, class T_VECTOR >
            void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );
          //! Finds scaling factor for that conf, collapsed dimension, and rank.
          typename t_Vector::value_type factor( size_t _kv, size_t _r, size_t _d );
          //! \brief Creates an X vector for fitting, for a single rank.
          //! \details Before the scaling is done, all ranks equal.
          template< class T_VECTOR >
            void create_X( size_t _i, T_VECTOR &_out );


          //! The configurations, arranged in columns.
          boost::shared_ptr<t_Configurations> configurations_;
          //! Holds current dimension being fitted.
          size_t dim;
          //! Pointer to separable function being minimized.
          Indirection::Pointer<t_Separables> separables_;
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
    std::ostream& operator<<( std::ostream& _stream, const Collapse<T_TRAITS>& _col )
      { return _stream << _col.separables() << "\n"; } 

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
