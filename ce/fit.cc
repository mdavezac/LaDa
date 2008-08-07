//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<iomanip>
#include<algorithm>
#include<numeric>
#include<iostream>
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/filesystem/path.hpp>
#include<boost/filesystem/operations.hpp>
#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>

#include <opt/random.h>
#include <opt/algorithms.h>
#include <opt/gsl_mins.h>

#include "fit.h"

namespace CE
{
  void Fit :: add_to_A_n_b( t_Vector &_A, t_Vector &_b,
                            const t_StructPis &_pis,
                            const types::t_real _weight,
                            const types::t_real _energy )
  {
    __ASSERT( i_pis->size() != nb_cls, "Inconsistent number of targets pis.\n" )
    // loop over alpha.
    t_ESums :: iterator i_b = _b.begin();
    t_PSums :: array_type :: iterator i_A = _A.data().begin();
    t_StructPis :: const_iterator i_pi_begin = _pis.begin();
    t_StructPis :: const_iterator i_pi_end = i_pis.end();
    t_StructPis :: const_iterator i_alphapi( i_pi_begin );
    for(; i_alphapi != i_pi_end; ++i_alphapi, ++i_esum)
    {
      __ASSERT( i_esum == _b.end(), "Iterator out of range.\n" )
      *i_b += _weight * energy * (*i_alphapi);

      // loop over betas.
      t_StructPis :: const_iterator i_betapi( i_pi_begin );
      for(; i_betapi != i_pi_end; ++i_betapi, ++i_A )
      {
        __ASSERT( i_A == _A.data().end(), "Iterator out of range.\n" )
        *i_A += _weight * (*i_alphapi) * (*i_betapi);
      }
    } // end of loop over alphas.
  } 

  void Fit :: create_A_n_b( t_Vector &_A, t_Vector &_b )
  {
    __DEBUGTRYBEGIN
    namespace bl = boost::lambda;
    __ASSERT( pis.size() != structures.size(),
              "Inconsistent number of structures and pis.\n" )
    __ASSERT( structures.size() != weights.size(),
              "Inconsistent number of structures and weights.\n" )

    // Resizes the clusters and fills with zeros.
    _b.resize( nb_cls );
    _A.resize( nb_cls, nb_cls );
    std::fill( _b.begin(), _b.end(), 0e0 );
    std::fill( _A.data().begin(), _A.data().end(), 0e0 );

    // loop over targets.
    t_Structures :: const_iterator i_target = structures.begin();
    t_Structures :: const_iterator i_target_end = structures.end();
    t_Pis :: const_iterator i_pis = pis.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    for( types::t_unsigned i(0); 
         i_target != i_target_end; 
         ++i_target, ++i_pis, ++i_w, ++i )
    {
      if( std::find( excluded.begin(), excluded.end(), i ) != excluded.end() )
        continue; 

      add_to_A_n_b( _A, _b, *i_pis, *i_w, i_target->energy );
    } // end of loop over betas.
    __DEBUGTRYEND(, "Error in Fit::create_A_n_b.\n" )
  } 

} // end of namespace CE

