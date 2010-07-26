#include "LaDaConfig.h"

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

#include <math/random.h>
#include <opt/algorithms.h>
#include <minimizer/interpolated_gradient.h>

#include "regularization.h"

namespace LaDa
{
  namespace CE
  {
    Regulated::t_Return Regulated :: operator()( const t_Vector& _arg ) const
    {
      LADA_DEBUG_TRY_BEGIN
      t_Vector x(nb_cls );
      std::fill( x.begin(), x.end(), 0e0 );
      std::pair< opt::ErrorTuple, opt::ErrorTuple >
        result = leave_one_out( *static_cast<const t_Base*>(this), cgs, x, &_arg[0], false );

      return result.second.variance();
      LADA_DEBUG_TRY_END(, "Error in Regulated::operator().\n" )
    }

    Regulated::t_Return Regulated :: operator()( const t_Arg& _arg ) const
    {
      LADA_DEBUG_TRY_BEGIN
      t_Vector x(nb_cls );
      std::fill( x.begin(), x.end(), 0e0 );
      std::pair< opt::ErrorTuple, opt::ErrorTuple >
        result = leave_one_out( *static_cast<const t_Base*>(this), cgs, x, &_arg[0], false );

      return result.second.variance();
      LADA_DEBUG_TRY_END(, "Error in Regulated::operator().\n" )
    }

    void Regulated :: gradient( const t_Arg& _arg,
                                types::t_real *_gradient ) const
    {
      LaDa::Minimizer::interpolated_gradient( *this, _arg, cgs, _gradient, 1, 1e-2 ); 
    }

    Cluster Regulated :: reduce()
    {
      LADA_DEBUG_TRY_BEGIN
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
        LADA_NASSERT( pis.size()-1 == clusters.size(),
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
      LADA_DEBUG_TRY_END(, "Error in Regulated::reduce().\n" )
    }

    opt::ErrorTuple Regulated :: fit( t_Arg& _arg, const types::t_real *_weights ) const
    {
      t_Vector arg( _arg.size() ); 
      std::copy( _arg.begin(), _arg.end(), arg.begin() );
      return Regulated::fit( arg, _weights );
    }
    std::pair< opt::ErrorTuple, opt::ErrorTuple > 
      Regulated :: loo( const types::t_real *_weights ) const
      {
        t_Vector x( nb_cls ); 
        std::fill( x.begin(), x.end(), 0 );
        return leave_one_out( *static_cast<const t_Base*>(this), cgs, x, _weights, false );
      }
  } // end of namespace CE
}
