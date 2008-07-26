//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<algorithm>
#include<boost/lambda/lambda.hpp>

#include "regularization.h"

namespace CE

  void Regulated :: init_sums()
  {
    namespace bl = boost::lambda;
    __ASSERT( pis.size() != structures.size(),
              "Inconsistent number of structures and pis.\n" )
    __ASSERT( structures.size() != weights.size(),
              "Inconsistent number of structures and weights.\n" )

    // Resizes the clusters and fills with zeros.
    size_t nb_cls( clusters.size() );
    numsums.resize( nb_cls );
    denumsums.resize( nb_cls );
    std::fill( numsums.begin(), numsums.end(), 0e0 );
    t_DenumSums :: iterator i_denums = denumsums.begin();
    t_DenumSums :: iterator i_denums_end = denumsums.begin();
    for(; i_denums != i_denums_end; ++i_denum )
    {
      i_denums->resize();
      std::fill( i_denums->begin(), i_denums_end(), 0e0 );
    }

    // loop over structures.
    t_Structures :: const_iterator i_str = structures.begin();
    t_Structures :: const_iterator i_str_end = structures.end();
    t_Pis :: const_iterator i_pis = pis.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    for(; i_str != i_str_end; ++i_str, ++i_pis, ++i_w )
    {
      __ASSERT( i_pis->size() == nb_cls, "Inconsistent number of structural pis.\n" )
      // loop over alpha.
      t_NumSums :: iterator i_num = numsums.begin();
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      i_denums = denumsums.begin();
      for(; i_denums != i_denums_end; ++i_num, ++i_denums, ++i_alphapi )
      {
        *i_num += (*i_w) * i_str->energy * (*i_alphapi);

        // loop over betas.
        t_DenumSums :: value_type :: iterator i_denum = i_denum_alpha->begin();
        t_DenumSums :: value_type :: iterator i_denum_end = i_denum_alpha->end();
        t_StructPis :: const_iterator i_betapi = i_pis->begin();
        for(; i_denum != i_denum_end; ++i_denum, ++i_beta )
          *i_denum += (*i_w) * (*i_alphapi) * (*i_betapi);

      } // end of loop over alphas.
    } // end of loop over betas.
  } 

  void operator( const t_Arg &_arg, t_Return &_return )
  {
    __ASSERT( _arg.size() != clusters.size(),
              "Inconsistent number of arguments and clusters.\n" )
    __ASSERT( structures.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( structures.size() != pis.size(),
              "Inconsistent number of returns and pis.\n" )

    if( _return.size() != structures.size() ) 
      _return.resize( structures.size() );

    // loop over return values.
    t_Return :: iterator i_ret = _return.begin();
    t_Return :: iterator i_ret_end = _return.end();
    t_Structures :: const_iterator i_str = structures.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for(; i_ret != i_ret_end; ++i_ret, ++i_str, ++i_w, ++i_pis)
    {
      __ASSERT( i_pis->size() == numsums.size(),
                "Inconsistent number of pis and numsums.\n" )
      __ASSERT( i_pis->size() == denumsums.size(),
                "Inconsistent number of pis and denumsums.\n" )
      *i_ret = (*i_w) * i_str->energy;

      // loop over alpha.
      t_Arg :: const_iterator i_arg = _args.begin();
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      t_DeNumSums :: const_iterator i_denums = denumsums.begin();
      for(; i_alphapi != i_pi_end; ++i_alphapi, ++i_denums )
      {
        // loop over beta.
        t_StructPis :: const_iterator i_betapi = i_pis->begin();
        t_DeNumSums :: const_iterator i_denums = denumsums.begin();
        t_NumSums :: const_iterator i_num = numsums.begin();
        for(; i_betapi != i_pi_end; ++i_betapi, ++i_num, ++i_denums )
        {
          types::t_real num, denum;
          num  = (*i_num) - (*i_w) * i_str->energy * (*i_betapi );
          num *= (*i_alphapi);
          denum = ( *i_denum - (*i_w) * (*i_alphapi) * (*i_betapi) );
          if( i_alphapi == i_betapi ) denum += (*i_arg) * (*i_arg);
          
          *i_ret += num / denum;
        } // end of loop over beta.
      }  // end of loop over alpha.
    } // end of loop over return values.
  }

  void jacobian( const t_Arg &_arg, t_Jacobian &_jacobian )
  {
    __ASSERT( _arg.size() != clusters.size(),
              "Inconsistent number of arguments and clusters.\n" )
    __ASSERT( structures.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( structures.size() != pis.size(),
              "Inconsistent number of returns and pis.\n" )

    if( _return.size() != structures.size() ) 
      _return.resize( structures.size() );

    // loop over return values.
    t_Return :: iterator i_ret = _return.begin();
    t_Return :: iterator i_ret_end = _return.end();
    t_Structures :: const_iterator i_str = structures.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for(; i_ret != i_ret_end; ++i_ret, ++i_str, ++i_w, ++i_pis)
    {
      __ASSERT( i_pis->size() == numsums.size(),
                "Inconsistent number of pis and numsums.\n" )
      __ASSERT( i_pis->size() == denumsums.size(),
                "Inconsistent number of pis and denumsums.\n" )
      *i_ret = (*i_w) * i_str->energy;

      // loop over alpha.
      t_Arg :: const_iterator i_arg = _args.begin();
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      t_DeNumSums :: const_iterator i_denums = denumsums.begin();
      for(; i_alphapi != i_pi_end; ++i_alphapi, ++i_denums )
      {
        // loop over beta.
        t_StructPis :: const_iterator i_betapi = i_pis->begin();
        t_DeNumSums :: const_iterator i_denums = denumsums.begin();
        t_NumSums :: const_iterator i_num = numsums.begin();
        for(; i_betapi != i_pi_end; ++i_betapi, ++i_num, ++i_denums )
        {
          types::t_real num, denum;
          num  = (*i_num) - (*i_w) * i_str->energy * (*i_betapi );
          num *= (*i_alphapi);
          denum = ( *i_denum - (*i_w) * (*i_alphapi) * (*i_betapi) );
          if( i_alphapi == i_betapi ) denum += (*i_arg) * (*i_arg);
          
          *i_ret += num / denum;
        } // end of loop over beta.
      }  // end of loop over alpha.
    } // end of loop over return values.
  }


} // end of namespace CE

