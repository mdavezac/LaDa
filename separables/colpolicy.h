#ifndef _CE_COLLAPSE_POLICY_H_
#define _CE_COLLAPSE_POLICY_H_

#include "LaDaConfig.h"

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/indirection.h>

#include "prepare.h"

namespace LaDa
{
  namespace CE
  {
    namespace Policy
    {
  //   //! \Initialization Policy for mixed approach.
  //   class Initialization
  //   {
  //      //! Functor.
  //      template< class T_THIS >
  //      static void init( const size_t &_rank, const size_t &_dimensions, 
  //                        T_THIS& _mixed,  );
  //   };
      
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
            LowMemUpdate() {}
            //! Destructor.
            virtual ~LowMemUpdate() {}
            
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
            void operator()( const t_Separables &_sep,
                             const t_Configurations &_confs );
            //! Updates only dimension \a _dim. 
            void operator()( size_t _dim,
                             const t_Separables &_sep,
                             const t_Configurations &_confs )
            //! Updates only dimension \a _dim. 
              { operator()( _sep, _confs ); }
            //! Returns the factor for configurations _kv, and rank _r, and dimension.
            typename t_Vector :: value_type
              factor( size_t _kv, size_t _r, size_t _dim,
                      const t_Separables &_sep,
                      const t_Configurations &_confs ) const;

          protected:
            //! Recomputes factor from nothing.
            typename t_Vector :: value_type
              factor_from_scratch( size_t _kv, size_t _r, size_t _dim,
                                   const t_Separables &_sep,
                                   const t_Configurations &_confs ) const;


            //! Holds the sep. func. split up by rank and confs.
            t_CMatrix scales_;
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
            HighMemUpdate() {}
            //! Updates all dimension.
            void operator()( const t_Separables &_sep,
                             const t_Configurations &_confs );
            //! Updates only dimension \a _dim. 
            void operator()( size_t _dim,
                             const t_Separables &_sep,
                             const t_Configurations &_confs );
            //! Returns the factor for configurations _kv, and rank _r, and dimension.
            typename t_Vector :: value_type
              factor( size_t _kv, size_t _r, size_t _dim,
                      const t_Separables &_sep,
                      const t_Configurations &_confs ) const;

          protected:
            //! Updates scales_ from dimsplit.
            void update_scales( const t_Separables &_sep );
            //! Returns the factor for configurations _kv, and rank _r, and dimension.
            typename t_Vector :: value_type
              factor_from_scratch( size_t _kv, size_t _r, size_t _dim ) const;

            //! Holds the sep. func. split along confs, ranks, and dimensions.
            std::vector< t_CMatrix > dimsplit_;
            using t_Base :: scales_;
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

          //! Would modify A matrix and b vector.
          template< class T_MATRIX, class T_VECTOR >
            void operator()( const t_Separables &_sep, T_MATRIX &, T_VECTOR &, size_t _dim ) const {}

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
          BadRegularization() : NoReg<T_SEPARABLES>(), lambda(0) {}
          //! Copy Constructor.
          BadRegularization   ( const BadRegularization &_c )
                            : NoReg<T_SEPARABLES>( _c ), lambda( _c.lambda ) {}

          //! Would modify A matrix and b vector.
          template< class T_MATRIX, class T_VECTOR >
            void operator()( const t_Separables& _sep, T_MATRIX &,
                             T_VECTOR &, size_t _dim ) const;
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
          Regularization() : NoReg<T_SEPARABLES>(), lambda(0) {}
          //! Copy Constructor.
          Regularization   ( const Regularization &_c )
                         : NoReg<T_SEPARABLES>( _c ), lambda( _c.lambda ) {}
          //! Copy Constructor.
          Regularization   ( const NoReg<T_SEPARABLES> &_c )
                         : NoReg<T_SEPARABLES>( _c ), lambda( 0 )  {}
          //! Assignement Operator.
          template< class TT_SEPARABLES >
            void operator=( const Regularization<TT_SEPARABLES> & _c ) { lambda = _c.lambda; }

          //! Would modify A matrix and b vector.
          template< class T_MATRIX, class T_VECTOR >
            void operator()( const t_Separables &_sep, T_MATRIX &, T_VECTOR &, size_t _dim ) const;
      };
    } // end of Policy namespace.

  }
} // namespace LaDa
#include "colpolicy.impl.h"

#endif
