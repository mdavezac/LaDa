//
//  Version: $Id$
//
#ifndef _CE_COLLAPSE_POLICY_H_
#define _CE_COLLAPSE_POLICY_H_

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
          //! Helps to obtain a differently typed update policy.
          template< class T_SEP, class T_MAP, class T_C = T_CONFS > struct rebind
          { 
            //! Type of the regularization policy.
            typedef LowMemUpdate< T_SEP, T_MAP, T_C > type;
          };
          //! Type of the separable functions.
          typedef T_SEPARABLES t_Separables;
          //! Type of the mapping from configurations to configurations.
          typedef T_MAPPING t_Mapping;
          //! Type of the mapping from configurations to configurations.
          typedef T_CONFS t_Configurations;
          //! Type of the vectors.
          typedef typename t_Separables :: t_Vector t_Vector;
          //! Type of the matrices.
          typedef typename t_Separables :: t_Matrix t_Matrix;
          //! Type of the constructible matrices.
          typedef boost::numeric::ublas::matrix
                    < typename t_Matrix :: value_type > t_CMatrix;

          //! the constructor.
          LowMemUpdate   ( const t_Mapping &_mapping )
                       : mapping_(_mapping) {}
          //! Destructor.
          ~LowMemUpdate() {}
          
          void print_scales() const
          {
            namespace bblas = boost::numeric::ublas;
            std::cout << "scales_ " << (long int) (&scales_) << "\n";
           //for( size_t i(0); i < configurations_->size2(); ++i )
           //{
           // bblas::matrix_column<const t_CMatrix> scaling( scales_, i );
           // std::cout << "3scaling " << i << ": " << scaling << "\n";
           //}
          }

          //! Updates all dimension.
          void operator()();
          //! Updates only dimension \a _dim. 
          void operator()( size_t _dim ) { operator()(); }
          //! Returns the factor for configurations _kv, and rank _r, and dimension.
          typename t_Vector :: value_type
            factor( size_t _kv, size_t _r, size_t _dim) const;
          //! Updates pointer to separable function.
          void init( const t_Separables& _sep );
          //! Updates pointer to configurations.
          void init( const boost::shared_ptr<t_Configurations> &_confs )
            { configurations_ = _confs; } 

        protected:
          //! Recomputes factor from nothing.
          typename t_Vector :: value_type
            factor_from_scratch( size_t _kv, size_t _r, size_t _dim) const;

          //! Holds the sep. func. split up by rank and confs.
          t_CMatrix scales_;
          //! A reference to the config. mapping.
          const t_Mapping &mapping_;
          //! A reference to the config. mapping.
          boost::shared_ptr<t_Configurations> configurations_;
          //! A reference to the separable function.
          const t_Separables * separables_;
      };

    template< class T_SEPARABLES, class T_MAPPING, class T_CONFS >
      class HighMemUpdate : protected LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS >
      {
        public:
          //! Helps to obtain a differently typed update policy.
          template< class T_SEP, class T_MAP, class T_C = T_CONFS > struct rebind
          { 
            //! Type of the regularization policy.
            typedef HighMemUpdate< T_SEP, T_MAP, T_C > type;
          };
          //! Type if the base class.
          typedef LowMemUpdate<T_SEPARABLES, T_MAPPING, T_CONFS > t_Base;
          //! Type of the separable functions.
          typedef T_SEPARABLES t_Separables;
          //! Type of the mapping from configurations to configurations.
          typedef T_MAPPING t_Mapping;
          //! Type of the mapping from configurations to configurations.
          typedef T_CONFS t_Configurations;
          //! Type of the vectors.
          typedef typename t_Separables :: t_Vector t_Vector;
          //! Type of the matrices.
          typedef typename t_Separables :: t_Matrix t_Matrix;
          //! Type of the constructible matrices.
          typedef boost::numeric::ublas::matrix
                    < typename t_Matrix :: value_type > t_CMatrix;

          //! Constructor.
          HighMemUpdate   ( const t_Mapping &_mapping )
                        : t_Base( _mapping ) {}
          //! Updates all dimension.
          void operator()();
          //! Updates only dimension \a _dim. 
          void operator()( size_t _dim );
          //! Returns the factor for configurations _kv, and rank _r, and dimension.
          typename t_Vector :: value_type
            factor( size_t _kv, size_t _r, size_t _dim) const;
          //! Updates pointer to separable function.
          void init( const t_Separables& _sep );
          //! Updates pointer to configurations.
          void init( const boost::shared_ptr<t_Configurations> &_confs )
            { t_Base::init( _confs ); }

        protected:
          //! Updates scales_ from dimsplit.
          void update_scales();
          //! Returns the factor for configurations _kv, and rank _r, and dimension.
          typename t_Vector :: value_type
            factor_from_scratch( size_t _kv, size_t _r, size_t _dim) const;

          //! Holds the sep. func. split along confs, ranks, and dimensions.
          std::vector< t_CMatrix > dimsplit_;
          using t_Base :: scales_;
          using t_Base :: mapping_;
          using t_Base :: configurations_;
          using t_Base :: separables_;
      };

    template< class T_SEPARABLES >
    class NoReg 
    {
      template< class TT_SEPARABLES > friend class NoReg;
      public:
        //! Helps to obtain a differently typed regularization.
        template< class T_SEP > struct rebind
        { 
          //! Type of the regularization policy.
          typedef NoReg< T_SEP > type;
        };
        //! Type of the separable function.
        typedef T_SEPARABLES t_Separables;
        //! Constructor
        NoReg() : separables_(NULL) {};
        //! Constructor
        template <class TT_SEPARABLES> 
          NoReg( const NoReg<TT_SEPARABLES> &_c ) : separables_( _c.separables_ ){};
        //! Updates pointer to separable function.
        void init( const t_Separables& _sep ) { separables_ = &_sep; }

        //! Would modify A matrix and b vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &, T_VECTOR &, size_t _dim ) const {}

      protected:
        //! A reference to the separable function.
        const t_Separables *separables_;
    };
    template< class T_SEPARABLES >
    class BadRegularization : public NoReg<T_SEPARABLES>
    {
      public:
        //! Helps to obtain a differently typed regularization.
        template< class T_SEP > struct rebind
        { 
          //! Type of the regularization policy.
          typedef BadRegularization< T_SEP > type;
        };
        //! Type of the separable function.
        typedef T_SEPARABLES t_Separables;

        //! Regulation factor.
        types::t_real lambda;

        //! Constructor
        BadRegularization() : NoReg<T_SEPARABLES>(), lambda(0) {};
        //! Copy Constructor.
        template <class TT_SEPARABLES> 
          BadRegularization   ( const BadRegularization<TT_SEPARABLES> &_c )
                            : NoReg<T_SEPARABLES>( _c ), lambda( _c.lambda ) {};
        //! Copy Constructor.
        template <class TT_SEPARABLES> 
          BadRegularization   ( const NoReg<TT_SEPARABLES> &_c )
                            : NoReg<T_SEPARABLES>( _c ), lambda( 0 )  {};

        //! Would modify A matrix and b vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &, T_VECTOR &, size_t _dim ) const;

        using NoReg<T_SEPARABLES> :: init;
      protected:
        using NoReg<T_SEPARABLES> :: separables_;
    };
    template< class T_SEPARABLES >
    class Regularization : public NoReg<T_SEPARABLES>
    {
      public:
        //! Helps to obtain a differently typed regularization.
        template< class T_SEP > struct rebind
        { 
          //! Type of the regularization policy.
          typedef Regularization< T_SEP > type;
        };
        //! Type of the separable function.
        typedef T_SEPARABLES t_Separables;

        //! Regulation factor.
        types::t_real lambda;

        //! Constructor
        Regularization() : NoReg<T_SEPARABLES>(), lambda(0) {};
        //! Copy Constructor.
        template <class TT_SEPARABLES> 
          Regularization   ( const Regularization<TT_SEPARABLES> &_c )
                         : NoReg<T_SEPARABLES>( _c ), lambda( _c.lambda ) {};
        //! Copy Constructor.
        template <class TT_SEPARABLES> 
          Regularization   ( const NoReg<TT_SEPARABLES> &_c )
                         : NoReg<T_SEPARABLES>( _c ), lambda( 0 )  {};

        //! Would modify A matrix and b vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &, T_VECTOR &, size_t _dim ) const;

        using NoReg<T_SEPARABLES> :: init;
      protected:
        using NoReg<T_SEPARABLES> :: separables_;
    };
  } // end of Policy namespace.

}

#include "colpolicy.impl.h"

#endif
