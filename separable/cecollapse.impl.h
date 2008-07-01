//
//  Version: $Id$
//

#include <opt/random.h>

namespace Separable
{
  template<class T_FUNCTION> template< class T_VECTORS, class T_VECTOR, class T_EWWIGHTS >
  void Collapse<T_FUNCTION> :: init( const T_VECTORS &_x, const T_VECTOR &_w,
                                     const T_EWWIGHTS &_eweights ) 
  {
    __DEBUGTRYBEGIN

    t_Base :: init( _x, _w );

    __ASSERT( _x.size() != _eweights.size(), "Inconsistent array-sizes on input.\n" )
    eweights.resize( _eweights.size() );
    typename T_EWEIGHTS :: const_iterator _i_eweights = _eweights.begin();
    t_eWeights :: iterator i_eweights = eweights.begin();
    t_eWeights :: iterator i_eweights_end = eweights.begin();
    for(; i_eweights != i_eweights_end; ++i_eweights, ++_i_eweights )
    {
      i_eweights->resize( _i_eweights->size() );
      std::copy( _i_eweights->begin(), _i_eweights->end(), i_eweights->begin() );
    }

    __DEBUGTRYEND(, "Error while initializing Collapse function for equivalent structures.\n" )
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
        __ASSERT( _targets.size() != nb_targets,
                  "Unexpected array-size on input.\n" )
        __ASSERT( weights.size() != nb_targets,
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
       
        if( ( not is_initialized ) or ( ( not do_update ) and _dim == 0 ) )
          initialize_factors( _coefs );
        else if ( do_update ) update_factors( _dim, _coefs );
        is_initialized = true;   
       
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
                     i_starget < i_starget_end;
                     ++o, ++i_iexpI, ++i_iexpII, ++i_starget, ++i_weight )
                {
                  __ASSERT( i_iexpI == i_rexpI->end(), "Iterator out of range.\n" )
                  __ASSERT( factors.size() <= rI * nb_targets + o, 
                            "Unexpected array-size.\n" )
                  __ASSERT( factors.size() <= rII * nb_targets + o, 
                            "Unexpected array-size.\n" )
       
                  // Computes first factor.
                  t_Type UI( 1 );
                  typename t_Factors :: value_type :: const_iterator
                     i_fac = factors[ rI * nb_targets + o].begin();
                  for(types::t_unsigned d(0); d < D; ++i_fac, ++d )
                    if( d != _dim )  UI *= (*i_fac);
                 
                  // Computes second factor.
                  t_Type UII( 1 ); i_fac = factors[ rII * nb_targets + o ].begin();
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


  template< class T_FUNCTION > template< class T_VECTORS, class T_VECTOR > 
    typename T_FUNCTION :: t_Return
      Collapse<T_FUNCTION> :: evaluate( const T_VECTORS &_coefs,
                                        const T_VECTOR &_target ) 
      {
#       ifdef _LADADEBUG
          try {
#       endif
        initialize_factors( _coefs );

        typename t_Factors :: const_iterator i_facs = factors.begin();
        std::vector< t_Type > values( _target.size(), 0 );
        for( size_t r(0); r < nb_ranks; ++r )
        {
          typename std::vector< t_Type > :: iterator i_val = values.begin();
          for( size_t o(0); o < nb_targets; ++o, ++i_facs, ++i_val )
            *i_val += std::accumulate
                      (
                        i_facs->begin(), i_facs->end(), 1e0,
                        std::multiplies<t_Type>()
                      );
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
        return result;
#     ifdef _LADADEBUG
        } __CATCHCODE(, "Error while assessing squares.\n" )
#     endif
      }

}
