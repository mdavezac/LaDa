//
//  Version: $Id$
//

#include <opt/random.h>

namespace LaDa
{
  namespace Separable
  {
    template<class T_FUNCTION> template< class T_VECTORS, class T_VECTOR >
    void Collapse<T_FUNCTION> :: init( const T_VECTORS &_x, const T_VECTOR &_w ) 
    {
      weights.resize( _w.size() );
      std::copy( _w.begin(), _w.end(), weights.begin() );

      // finds maximum size of separable function.
      typename t_Function::t_Basis::const_iterator i_first = function.basis.begin();
      typename t_Function::t_Basis::const_iterator i_last = function.basis.end();
      for(D=0; i_first != i_last; ++i_first )
        if( D < i_first->basis.size() ) D = i_first->basis.size();
      // finds nb of target values
      nb_targets = _x.size();
      // finds nb of ranks
      nb_ranks = function.basis.size();

      // computes sizes
      sizes.resize( D );
      for( size_t d = 0; d < D; ++d )
      {
        typedef typename t_Function::t_Basis::const_iterator t_const_iterator;
        t_const_iterator i_rank = function.basis.begin();
        t_const_iterator i_rank_end = function.basis.end();
        for(; i_rank != i_rank_end; ++i_rank )
        {
          if( i_rank->basis.size() < d ) continue;
          sizes[d].push_back( i_rank->basis[d].basis.size() );
        }
      }

      // Computes expanded
      expanded.resize( D );
      typename t_Expanded :: iterator i_dexp = expanded.begin();
      t_Sizes :: const_iterator i_dsize = sizes.begin();
      for( types::t_unsigned d(0); d < D; ++d, ++i_dexp, ++i_dsize )
      {
        __ASSERT( i_dexp  == expanded.end(), "Iterator out of range.\n" )
        __ASSERT( i_dsize == sizes.end(),    "Iterator out of range.\n" )
        
        // prepares loop over ranks
        i_dexp->resize( i_dsize->size() );
        typename t_Expanded :: value_type :: iterator i_rexp = i_dexp->begin();
        typename t_Sizes :: value_type :: const_iterator i_rsize = i_dsize->begin();
        typename t_Sizes :: value_type :: const_iterator i_rsize_end = i_dsize->end();
        for( types::t_unsigned r(0); i_rsize != i_rsize_end;
             ++r, ++i_rexp, ++i_rsize ) 
        {
          __ASSERT( _x.size() != nb_targets, "Unexpected array size.\n" )
          __ASSERT( i_rexp  == i_dexp->end(),  "Iterator out of range.\n" )
          __ASSERT( i_rsize == i_dsize->end(), "Iterator out of range.\n" )
          __ASSERT( function.basis[r].basis[d].basis.size() != *i_rsize,
                    "Inconsistent size.\n" )

          // prepares loop over functions.
          i_rexp->resize( (*i_rsize) * nb_targets );
          typename t_Expanded :: value_type
                              :: value_type :: iterator i_iexp = i_rexp->begin();
          typename t_Function :: t_Basis :: value_type
                              :: t_Basis :: value_type 
                              :: t_Basis
                              :: const_iterator i_ifunc = function.basis[r]
                                                                  .basis[d]
                                                                  .basis.begin();
          // loop over basis functions.
          for( types::t_unsigned i(0); i < *i_rsize; ++i, ++i_ifunc )
            // loop over structural target values.
            for( types::t_unsigned o(0); o < nb_targets; ++o, ++i_iexp )
            {
              __ASSERT( _x[o].size() != D, "Unexpected array-size.\n" )
              __ASSERT( i_iexp == i_rexp->end(), "Iterator out of range.\n" )
              *i_iexp = (*i_ifunc)( _x[o][d] );
            }

          __ASSERT( i_iexp  != i_rexp->end(), "Unexpected iterator position.\n" )
          __ASSERT( i_ifunc != function.basis[r].basis[d].basis.end(),
                    "Unexpected iterator position.\n" )
        } // end of loop over ranks
        __ASSERT( i_rexp  != i_dexp->end(),  "Unexpected iterator position.\n" )
        __ASSERT( i_rsize != i_dsize->end(), "Unexpected iterator position.\n" )

      } // end of loop over dimensions.
      __ASSERT( i_dexp  != expanded.end(), "Unexpected iterator position.\n" )
      __ASSERT( i_dsize != sizes.end(),    "Unexpected iterator position.\n" )

    } // end of init member function.

    template<class T_FUNCTION>  
      template< class T_VECTOR, class T_MATRIX, class T_VECTORS >
        void Collapse<T_FUNCTION>::operator()( T_VECTOR &_b, T_MATRIX &_A, 
                                               types::t_unsigned _dim,
                                               const T_VECTOR &_targets,
                                               const T_VECTORS &_coefs     )
        {
          __DEBUGTRYBEGIN

          __ASSERT( _dim  >= D, "Input dimension out of range.\n" )
          __ASSERT( _coefs.size() != D,
                    "Inconsistent number of dimensions/coefficients.\n" )
          __ASSERT( _targets.size() != weights.size(),
                    "Unexpected array-size on input.\n" )
          __ASSERT( function.coefs.size() != nb_ranks,
                    "Unexpected array-size on input.\n" )
          __ASSERT
          ( 
            _coefs[_dim].size() != std::accumulate( sizes[_dim].begin(),
                                                    sizes[_dim].end(), 0 ),
            "Inconsistent number of ranks/dimensions/coefficients.\n" 
          )
          
          __ASSERT( sizes.size()  != D, "Inconsistent sizes.\n" )
          __ASSERT( factors.size() == _coefs[_dim].size(), "Inconsistent sizes.\n" )
          __ASSERT( expanded.size() <= _dim, "Inconsistent sizes.\n" )
       
          if( _b.size() != _coefs[_dim].size() )
          {
            _b.resize( _coefs[_dim].size() );
            _A.resize( _b.size(), _b.size() );
          }
          __ASSERT( _A.size1() != _b.size(), "Wrong size for matrix _A.\n" );
          __ASSERT( _A.size2() != _b.size(), "Wrong size for matrix _A.\n" );
       
          typename T_VECTOR :: iterator i_ctarget = _b.begin();
          typename t_Expanded :: value_type 
                              :: const_iterator i_rexpI = expanded[_dim].begin();
          typename t_Expanded :: value_type
                              :: const_iterator i_rexp_end = expanded[_dim].end();
          typename t_Sizes :: value_type
                           :: const_iterator i_sizeI = sizes[_dim].begin();
          typename t_Sizes :: value_type :: const_iterator i_size_first( i_sizeI );
          size_t k1(0);
   
          // loops over ranks I
          for( size_t rI(0); i_rexpI < i_rexp_end; ++rI, ++i_rexpI, ++i_sizeI ) 
          {
            // loop over basis functions I
            typename t_Expanded :: value_type :: value_type
                                :: const_iterator i_iexp_inc = i_rexpI->begin();
            for( size_t iI(0); iI < *i_sizeI;
                 ++iI, i_iexp_inc += nb_targets, ++i_ctarget, ++k1 )
            {
              __ASSERT( i_ctarget == _b.end(),        "Iterator out of range.\n" )
              __ASSERT( i_iexp_inc == i_rexpI->end(), "Iterator out of range.\n" )
              *i_ctarget = typename T_VECTOR :: value_type(0);
       
              // loops over ranks I
              typename t_Expanded :: value_type 
                                  :: const_iterator i_rexpII = expanded[_dim].begin();
              typename t_Sizes :: value_type 
                               :: const_iterator i_sizeII( i_size_first );
              for( size_t rII(0), k2(0); i_rexpII < i_rexp_end;
                   ++rII, ++i_rexpII, ++i_sizeII ) 
              {
                __ASSERT( i_sizeII == sizes[_dim].end(), "Iterator out of range.\n" )
       
                // loop over basis functions II
                typename t_Expanded :: value_type :: value_type
                                    :: const_iterator i_iexpII = i_rexpII->begin();
                for(size_t iII(0); iII < *i_sizeII; ++iII, ++k2) //, ++i_A)
                {
                  // Initializes element to 0 before summing over structural targets.
                  _A( k1, k2 ) = 0;
       
                  // loop over target values
                  typename t_Expanded :: value_type :: value_type
                                      :: const_iterator i_iexpI( i_iexp_inc );
                  typename T_VECTOR :: const_iterator i_starget = _targets.begin();
                  typename T_VECTOR :: const_iterator i_starget_end = _targets.end();
                  typename t_Weights :: const_iterator i_weight = weights.begin();
                  for( types::t_unsigned o(0);
                       i_starget != i_starget_end;
                       ++o, ++i_iexpI, ++i_iexpII, ++i_starget, ++i_weight )
                  {
                    __ASSERT( i_iexpI == i_rexpI->end(), "Iterator out of range.\n" )
                    __ASSERT( factors.size() <= rI * nb_targets + o, 
                              "Unexpected array-size.\n" )
                    __ASSERT( factors.size() <= rII * nb_targets + o, 
                              "Unexpected array-size.\n" )
       
                    // Computes first factor.
                    t_Type UI( function.coefs[ rI ] );
                    typename t_Factors :: value_type :: const_iterator
                       i_fac = factors[ rI * nb_targets + o].begin();
                    for(types::t_unsigned d(0); d < D; ++i_fac, ++d )
                      if( d != _dim )  UI *= (*i_fac);
                   
                    // Computes second factor.
                    t_Type UII( function.coefs[ rII ] ); 
                    i_fac = factors[ rII * nb_targets + o ].begin();
                    for(types::t_unsigned d(0); d < D; ++i_fac, ++d )
                      if( d != _dim ) UII *= (*i_fac);
       
                    // Computes matrix element.
                    _A(k1,k2) += (*i_iexpI) * UI * UII * (*i_iexpII ) * (*i_weight);
       
                    // bypasses collapsed target on second dimension of A.
                    if( i_sizeII != i_size_first or iII != 0 ) continue;
                   
                    // Computes collapsed target element.
                    __ASSERT( i_ctarget == _b.end(), "Iterator out of range.\n" )
                    *i_ctarget += (*i_iexpI) * UI * (*i_starget) * (*i_weight);
                  } // loop over structural target values.
       
                } // loop over basis functions II
       
              } // end of loop over ranks II
       
            } // loop over basis functions I
            __ASSERT( i_iexp_inc != i_rexpI->end(),
                      "Unexpected iterator position.\n" )
       
          } // end of loop over ranks I

          __DEBUGTRYEND(, "Error while creating collapsed matrix and target.\n" )
       
        } // end of functor.

    template< class T_FUNCTION> template< class T_VECTORS >
      void Collapse<T_FUNCTION>::update_all( T_VECTORS &_coefs )
      {
        namespace bl = boost::lambda;
        __DEBUGTRYBEGIN

        __ASSERT( nb_ranks != function.basis.size(),
                  "Inconsistent rank size.\n" )

        if( factors.size() != nb_ranks * nb_targets )
          factors.resize( nb_ranks * nb_targets );
        typename t_Factors :: iterator i_facs = factors.begin();
        for( types::t_unsigned r(0); r < nb_ranks; ++r )
        {
          // loop over target values.
          for(types::t_unsigned o(0); o < nb_targets; ++o, ++i_facs )
          {
            if( i_facs->size() != D ) i_facs->resize( D );
            typename t_Factors :: value_type :: iterator i_fac = i_facs->begin();
            for( types::t_unsigned d(0); d < D; ++d, ++i_fac )
            {
              __ASSERT( function.basis.size() <= r,
                           "Inconsistent rank size.: " 
                        << function.basis.size() << " " << r << ".n" )

              // if there less than d dimensions for this rank, the factor is
              // evaluated to 1, such that it will have no influence.
              if( function.basis[r].basis.size() <= d )
                { *i_fac = t_Type(1); continue; }

              types::t_unsigned acc = std::accumulate( sizes[d].begin(),
                                                       sizes[d].begin() + r, 0 );
              __ASSERT( _coefs.size()  <= d, "Inconsistent dimension size.\n" )
              __ASSERT( _coefs[d].size() < acc + sizes[d][r],
                           "Inconsistent size: " << _coefs[d].size()
                         << " < " << acc + sizes[d][r] << ".\n" )
              __ASSERT( expanded[d][r].size() !=  sizes[d][r] * nb_targets,
                           "Unexpected array size. " << expanded[d][r].size()
                        << " != " << sizes[d][r] * nb_targets << ".\n" )

              typedef typename t_Expanded :: value_type
                                          :: value_type :: const_iterator t_itexp;
              typedef typename T_VECTORS :: value_type :: iterator t_itcoefs;
              t_itcoefs i_coef = _coefs[d].begin() + acc;
              t_itexp i_iexp = expanded[d][r].begin() + o;

              // normalize the coefficients of function family in (r,d)
              types::t_real norm(0);
              std::for_each
              (
                i_coef, i_coef + sizes[d][r],
                bl::var(norm) += bl::_1 * bl::_1
              );
              if( Fuzzy::neq( norm, 0e0 ) )
              { 
                norm = std::sqrt( norm );
                function.coefs[r] *= norm; 
                norm = 1e0 / norm;
              }
              else norm = 1e0;

              // updates the factors.
              *i_fac = t_Type(0);
              for( types::t_unsigned i( sizes[d][r] ); i > 0;
                   --i, i_iexp += nb_targets, ++i_coef )
              {
                __ASSERT( i_iexp == expanded[d][r].end(), "Iterator out of range.\n" )
                __ASSERT( i_coef == _coefs[d].end(), "Iterator out of range.\n" )
                *i_coef *= norm;
                *i_fac += (*i_coef) * (*i_iexp ); 
              }
              
            } // end of loop over dimensions
            __ASSERT( i_fac != i_facs->end(), "Unexpected iterator position.\n" )
          } // end of loop over target values.
        } // end of loop over ranks
        __ASSERT( i_facs != factors.end(), "Unexpected iterator position.\n" )

        __DEBUGTRYEND(, "Error while creating factors.\n" )
      }

    template< class T_FUNCTION > template< class T_VECTORS >
      void Collapse<T_FUNCTION>::update( types::t_unsigned _dim, T_VECTORS &_coefs )
      {
        namespace bl = boost::lambda;
        __DEBUGTRYBEGIN

        __ASSERT( _dim >= D, "Dimension out of range on input.\n" ) 
        __ASSERT( sizes.size() != D, "Unexpected array-size.\n" ) 
        __ASSERT( nb_ranks != function.basis.size(),
                  "Inconsistent rank size.\n" )
        __ASSERT( factors.size() != nb_ranks * nb_targets,
                  "Inconsistent rank size.\n" )

        typename t_Factors :: iterator i_facs = factors.begin();
        for( types::t_unsigned r(0); r < nb_ranks; ++r )
        {
          __ASSERT( sizes[_dim].size() != nb_ranks, "Unexpected array-size.\n" ) 

          // loop over target values.
          for(types::t_unsigned o(0); o < nb_targets; ++o, ++i_facs )
          {
            __ASSERT( i_facs == factors.end(), "Iterator out of range.\n" )
            __ASSERT( i_facs->size() != D, "Unexpected array size.\n" )
            __ASSERT( function.basis.size() <= r,
                       "Inconsistent rank size.: " 
                      << function.basis.size() << " " << r << ".n" )
             
            // if there less than d dimensions for this rank, the factor is
            // evaluated to 1, such that it will have no influence.
            if( function.basis[r].basis.size() <= _dim )
              { (*i_facs)[_dim] = t_Type(1); continue; }
            
            types::t_unsigned acc = std::accumulate( sizes[_dim].begin(),
                                                     sizes[_dim].begin() + r, 0 );
            __ASSERT( _coefs.size()  <= _dim, "Inconsistent dimension size.\n" )
            __ASSERT( _coefs[_dim].size() < acc + sizes[_dim][r],
                      "Inconsistent size.\n" )
            __ASSERT( expanded[_dim][r].size() !=  sizes[_dim][r] * nb_targets,
                         "Unexpected array size. " << expanded[_dim][r].size()
                      << " != " << sizes[_dim][r] * nb_targets << ".\n" )

            typedef typename t_Expanded :: value_type
                                        :: value_type :: const_iterator t_itexp;
            typedef typename T_VECTORS :: value_type :: iterator t_itcoefs;
       
            t_itcoefs i_coef = _coefs[_dim].begin() + acc;
            t_itexp i_iexp = expanded[_dim][r].begin() + o;
            
            // normalize the coefficients of function family in (r,d)
            types::t_real norm(0);
            std::for_each
            (
              i_coef, i_coef + sizes[_dim][r],
              bl::var(norm) += bl::_1 * bl::_1
            );
            if( Fuzzy::neq( norm, 0e0 ) )
            { 
              norm = std::sqrt( norm );
              function.coefs[r] *= norm; 
              norm = 1e0 / norm;
            }
            else norm = 1e0;

            (*i_facs)[_dim] = t_Type(0);
            for( types::t_unsigned i( sizes[_dim][r] ); i > 0;
                 --i, i_iexp += nb_targets, ++i_coef )
            {
              __ASSERT( i_iexp == expanded[_dim][r].end(), 
                        "Iterator out of range.\n" )
              __ASSERT( i_coef == _coefs[_dim].end(), "Iterator out of range.\n" )
              *i_coef *= norm;
              (*i_facs)[_dim] += (*i_coef) * (*i_iexp ); 
            }
            
          } // end of loop over target values.
        } // end of loop over ranks
        __ASSERT( i_facs != factors.end(), "Inconsistent size.\n" )

        __DEBUGTRYEND(, "Error while updating factors.\n" )
      }

    template< class T_FUNCTION > template< class T_VECTORS >
      void Collapse<T_FUNCTION> :: reassign( const T_VECTORS& _solution ) const
      {
        __DEBUGTRYBEGIN

        typedef T_VECTORS t_Vectors;
        typedef typename t_Vectors :: const_iterator const_solit;
        typedef typename t_Vectors :: value_type :: const_iterator const_soldimit;
        const_solit i_dim = _solution.begin();
        const_solit i_dim_end = _solution.end();
        for(size_t d=0; i_dim != i_dim_end; ++i_dim, ++d ) // loop over dimensions.
        {
          const_soldimit i_coef = i_dim->begin();
          const_soldimit i_coef_ = i_coef;
          const_soldimit i_coef_end = i_dim->end();
          for(size_t r=0; i_coef != i_coef_end; ++r )
          {
            typedef typename t_Function :: t_Basis :: value_type
                                        :: t_Basis :: value_type t_factor_func;
            t_factor_func &facfunc = function.basis[r].basis[d];
            i_coef_ += facfunc.coefs.size();
            __ASSERT( i_coef_end - i_coef_ < 0,
                      "Inconsistent number of coefficients.\n" )
            std::copy( i_coef, i_coef_, facfunc.coefs.begin() );
            i_coef = i_coef_;
          }
        }

        __DEBUGTRYEND(, "Error while assigning solution coefs to function.\n" )
      }

    template< class T_FUNCTION > template< class T_VECTORS >
      void Collapse<T_FUNCTION> :: create_coefs( T_VECTORS& _coefs,
                                                 t_Type _howrandom ) const
      {
        // Initializes norms to 1
        std::fill( function.coefs.begin(), function.coefs.end(), t_Type(1) );

        __DEBUGTRYBEGIN
        typedef T_VECTORS t_Vectors;
        typedef typename t_Vectors :: iterator t_solit;
        _coefs.resize( sizes.size() );
        t_solit i_dim = _coefs.begin();
        t_solit i_dim_end = _coefs.end();
        t_Sizes :: const_iterator i_sizes = sizes.begin();
        types::t_real ndim( sizes.size() );

        // loop over dimensions
        for(; i_dim != i_dim_end; ++i_dim, ++i_sizes)
        {
          types::t_int ri = std::accumulate( i_sizes->begin(), i_sizes->end(), 0 );
          i_dim->resize( ri );
          t_Sizes :: value_type :: const_iterator i_size = i_sizes->begin();
          t_Sizes :: value_type :: const_iterator i_size_end = i_sizes->end();
          typename t_Vectors :: value_type :: iterator i_coef = i_dim->begin();
          // loop over ranks.
          for(size_t r(0); i_size != i_size_end; ++i_size, ++r )
          {
            // First creates a normalized random vector.
            namespace bl = boost::lambda;
            typename t_Vectors :: value_type :: iterator i_c( i_coef ); 
            t_Type norm;
            do
            {
              norm = t_Type(0);
              const t_Type range(_howrandom);
              std::for_each 
              (
                i_c, i_c + *i_size, 
                bl::var( norm ) += (
                                     bl::_1 =   bl::constant(range) 
                                              * bl::bind( &opt::random::rng )
                                              + bl::constant( t_Type(1e0 - 0.5*range) ),
                                     bl::_1 * bl::_1
                                   )
              );
            }
            while( Fuzzy::eq( norm, t_Type(0) ) );
            norm = std::sqrt(norm);
            function[r] = norm;
            norm = t_Type(1) / norm;
            std::for_each 
            (
              i_c, i_c + *i_size, 
              bl::_1 *= bl::constant(norm)
            );
            i_coef += *i_size;
            // Then adds random norm to coefficient
          } // end of loop over ranks.
        } // end of loop over dimensions.
        __DEBUGTRYEND(, "Error while assigning solution coefs to function.\n" )
      }

    template< class T_FUNCTION > template< class T_VECTORS, class T_VECTOR > 
      typename T_FUNCTION :: t_Return
        Collapse<T_FUNCTION> :: evaluate( const T_VECTORS &_coefs,
                                          const T_VECTOR &_target ) 
        {
          __DEBUGTRYBEGIN

          typename t_Factors :: const_iterator i_facs = factors.begin();
          std::vector< t_Type > values( _target.size(), 0 );
          typename t_Function :: t_Coefs :: const_iterator i_coef = function.coefs.begin();
          for( size_t r(0); r < nb_ranks; ++r, ++i_coef )
          {
            typename std::vector< t_Type > :: iterator i_val = values.begin();
            for( size_t o(0); o < nb_targets; ++o, ++i_facs, ++i_val )
              *i_val += std::accumulate
                        (
                          i_facs->begin(), i_facs->end(), 1e0,
                          std::multiplies<t_Type>()
                        ) * (*i_coef);
          } // end of loop over ranks.

          t_Type result(0);
          typename std::vector< t_Type > :: iterator i_val = values.begin();
          typename T_VECTOR :: const_iterator i_target = _target.begin();
          typename T_VECTOR :: const_iterator i_target_end = _target.end();
          typename t_Weights :: const_iterator i_weight = weights.begin();
          for(; i_target != i_target_end; ++i_target, ++i_val, ++i_weight )
          {
            t_Type intermed = *i_target - *i_val;
            result += intermed * intermed * (*i_weight);
          }
          return result / (types::t_real) _target.size();
        
          __DEBUGTRYEND(, "Error while assessing squares.\n" )
        }

    template< class T_FUNCTION >
      void Collapse<T_FUNCTION> :: reset() 
      {
        expanded.clear();
        factors.clear();
        sizes.clear();
      }
  }
} // namespace LaDa
