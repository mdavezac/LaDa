//
//  Version: $Id$
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

#include "regularization.h"

namespace LaDa
{
  namespace CE
  {
    Regulated::t_Return Regulated :: operator()( const types::t_real *_arg ) const
    {
      __DEBUGTRYBEGIN
      t_Vector x(nb_cls );
      std::fill( x.begin(), x.end(), 0e0 );
      std::pair< opt::ErrorTuple, opt::ErrorTuple >
        result = leave_one_out( *static_cast<const t_Base*>(this), cgs, x, _arg, false );

      return result.second.variance();
      __DEBUGTRYEND(, "Error in Regulated::operator().\n" )
    }

    void Regulated :: gradient( const types::t_real * _arg,
                                types::t_real *_gradient ) const
    {
      __DEBUGTRYBEGIN
      types::t_real delta = 1e-2; // * cgs.tolerance;
      types::t_real delta_inv = 2e0 / delta;
      types::t_real right, left, keep;
      const types::t_real *i_arg( _arg );
      std::vector< types::t_real > arg( nb_cls );
      std::copy( i_arg, i_arg + nb_cls, arg.begin() );
      for( size_t i(0); i < nb_cls; ++i )
      {
        types::t_real keep( arg[i] );
        arg[i] += delta; 
        types::t_real left( Regulated::operator()( &arg[0] ) ); 
        arg[i] = keep - delta; 
        types::t_real right( Regulated::operator()( &arg[0] ) ); 
        *( _gradient + i ) = ( left - right ) * delta_inv;
      }
      __DEBUGTRYEND(, "Error in Regulated::gradient().\n" )
    }

    Cluster Regulated :: reduce()
    {
      __DEBUGTRYBEGIN
      namespace bblas = boost::numeric::ublas;
      t_Clusters :: iterator i_cls = clusters.begin(); 
      t_Clusters :: iterator i_cls_end = clusters.end(); 
      t_Clusters :: iterator i_found( i_cls );
      types::t_real least( -1e0 );
      types::t_unsigned index(0); ++i_cls;
      for(types::t_unsigned i(0); i_cls != i_cls_end; ++i_cls, ++i )
      {
        if( i_cls->front().size() < 2 ) continue;
        if( least > std::abs( i_cls->front().eci ) or least < 0e0 )
        {
          i_found = i_cls;
          least = std::abs(i_cls->front().eci );
          index = i;
        }
      }

      Cluster result( i_found->front() );
      clusters.erase( i_found );
      --nb_cls;

      std::vector<t_StructPis> :: iterator i_pis = pis.begin();
      std::vector<t_StructPis> :: iterator i_pis_end = pis.end();
      for(; i_pis != i_pis_end; ++i_pis )
      {
        __ASSERT( pis.size()-1 == clusters.size(),
                  "Inconsistent sizes of Pis and clusters.\n")

        t_Vector temp( *i_pis );
        i_pis->resize( nb_cls );
        t_Vector :: const_iterator i_temp = temp.begin();
        t_Vector :: const_iterator i_temp_end = temp.end();
        t_Vector :: const_iterator i_pi = i_pis->begin();
        for( size_t i(0); i_temp != i_temp_end; ++i_temp, ++i )
        {
          if( i == index ) continue;
          *i_pi == *i_temp;
          ++i_pi;
        }
      }

      return result;
      __DEBUGTRYEND(, "Error in Regulated::reduce().\n" )
    }
  } // end of namespace CE
}
