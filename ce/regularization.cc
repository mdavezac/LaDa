//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<iomanip>
#include<algorithm>
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>

#include <opt/random.h>
#include <opt/algorithms.h>
#include <opt/gsl_mins.h>

#include "regularization.h"

namespace CE
{
  namespace details
  {
    types::t_real innerprod( const boost::numeric::ublas::vector<types::t_real> &_a1,
                             const boost::numeric::ublas::vector<types::t_real> &_a2 )
      { return boost::numeric::ublas::inner_prod( _a1, _a2 ); }
  }
  void Regulated :: init_sums()
  {
    namespace bl = boost::lambda;
    __ASSERT( pis.size() != targets.size(),
              "Inconsistent number of targets and pis.\n" )
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of targets and weights.\n" )

    // Resizes the clusters and fills with zeros.
    esums.resize( nb_cls );
    psums.resize( nb_cls, nb_cls );
    std::fill( esums.begin(), esums.end(), 0e0 );
    std::fill( psums.data().begin(), psums.data().end(), 0e0 );

    // loop over targets.
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.end();
    t_Pis :: const_iterator i_pis = pis.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    for(; i_target != i_target_end; ++i_target, ++i_pis, ++i_w )
    {
      __ASSERT( i_pis->size() != nb_cls, "Inconsistent number of targets pis.\n" )
      // loop over alpha.
      t_ESums :: iterator i_esum = esums.begin();
      t_PSums :: array_type :: iterator i_psum = psums.data().begin();
      t_StructPis :: const_iterator i_pi_begin = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      t_StructPis :: const_iterator i_alphapi( i_pi_begin );
      for(; i_alphapi != i_pi_end; ++i_alphapi, ++i_esum )
      {
        __ASSERT( i_esum == esums.end(), "Iterator out of range.\n" )
        *i_esum += (*i_w) * (*i_target) * (*i_alphapi);

        // loop over betas.
        t_StructPis :: const_iterator i_betapi( i_pi_begin );
        for(; i_betapi != i_pi_end; ++i_betapi, ++i_psum )
        {
          __ASSERT( i_psum == psums.data().end(), "Iterator out of range.\n" )
          *i_psum += (*i_w) * (*i_alphapi) * (*i_betapi);
        }
      } // end of loop over alphas.
    } // end of loop over betas.
  } 

  Regulated::t_Return Regulated :: operator()( const types::t_real * _arg ) const
  {
    namespace bblas = boost::numeric::ublas;
    namespace bl = boost::lambda;
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( targets.size() != fittingpairs.size(),
              "Incorrect number of fitting pairs.\n" )

    types::t_real result(0), dummy;
//   opt::concurrent_loop
//   (
//     _arg, _arg + nb_cls, targets.begin(), 
//     bl::var(result) += 
//     (
//       bl::var(dummy) = bl::_1 - bl::_2,
//       bl::var(dummy) * bl::var(dummy) 
//     )
//   );
//   return result;

    t_Matrix A(nb_cls, nb_cls); 
    t_Vector x(nb_cls );
    // loop over target values.
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.end();
    t_Targets :: const_iterator i_w = weights.begin();
    t_FittingPairs :: const_iterator i_pair = fittingpairs.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for( ; i_target != i_target_end;
         ++i_target, ++i_w, ++i_pair, ++i_pis )
    {
      __ASSERT( i_w == weights.end(), "Iterator out of range.\n" )
      __ASSERT( i_pair == fittingpairs.end(), "Iterator out of range.\n" )
      __ASSERT( i_pis == pis.end(), "Iterator out of range.\n" )

      // construct matrix with weights.
      A = i_pair->first;
      const types::t_real *i_arg = _arg;
      for( size_t i(0); i < nb_cls; ++i, ++i_arg )
      {
        A(i,i) += (*i_arg) * (*i_arg);
      }
      const types::t_real range(0);
      std::for_each
      ( 
        x.begin(), x.end(),
        bl::_1 =  bl::bind( &opt::random::rng ) * range - range * 0.5e0
      );

      // fits.
      cgs( A, x, i_pair->second );

      // computes error for this target.
      types::t_real intermed( *i_target - bblas::inner_prod( x, *i_pis ) );
      // adds to result
      result += (*i_w) * intermed * intermed;
    } // end of loop over target values.

    return result / types::t_real( targets.size() );
  }

  void Regulated :: gradient( const types::t_real * _arg,
                              types::t_real *_gradient ) const
  {
    types::t_real delta = 1e-4; // * cgs.tolerance;
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
  }

  void Regulated :: init( const t_Structures &_structures )
  {
    namespace bl = boost::lambda;
    // inits weights
    weights.resize( _structures.size() );
    std::fill( weights.begin(), weights.end(), 1e0 );

    // inits targets.
    targets.resize( _structures.size() );
    std::transform
    ( 
      _structures.begin(), _structures.end(), targets.begin(),
      bl::bind( &t_Structures::value_type::energy, bl::_1 )
    );

    // init pis.
    find_pis( clusters, _structures, pis );
    // t_Structures :: const_iterator i_structure = _structures.begin();
    // t_Structures :: const_iterator i_structure_end = _structures.end();
    // t_Pis ::const_iterator i_pis = pis.begin();
    // for(; i_structure != i_structure_end; ++i_structure, ++i_pis )
    // {
    //   std::cout << i_structure->name << ": ";
    //   std::for_each
    //   (
    //     i_pis->begin(), i_pis->end(),
    //     bl::var( std::cout ) << std::fixed << bl::_1 << " " 
    //   );
    //   std::cout << "\n";
    // }

    // initializes intermediate quantities.
    nb_cls = clusters.size();
    init_sums();
    construct_pairs();
  } 
 
  types::t_unsigned Regulated :: reduce()
  {
    namespace bblas = boost::numeric::ublas;
    t_Clusters :: iterator i_cls = clusters.begin(); 
    t_Clusters :: iterator i_cls_end = clusters.end(); 
    t_Clusters :: iterator i_found( i_cls );
    types::t_real least( std::abs(i_cls->front().eci ) );
    types::t_unsigned index(0); ++i_cls;
    for(types::t_unsigned i(0); i_cls != i_cls_end; ++i_cls, ++i )
      if( least > std::abs(i_cls->front().eci ) )
      {
        i_found = i_cls;
        least = std::abs(i_cls->front().eci );
        index = i;
      }

    clusters.erase( i_found );
    --nb_cls;

    t_Pis :: iterator i_pis = pis.begin();
    t_Pis :: iterator i_pis_end = pis.end();
    for(; i_pis != i_pis_end; ++i_pis )
    {
      __ASSERT( pis.size()-1 == clusters.size(),
                "Inconsistent sizes of Pis and clusters.\n")

      t_Vector temp( nb_cls );
      if( index != 0 ) 
        bblas::noalias( bblas::subrange( temp, 0, index ) )
           = bblas::subrange( *i_pis, 0, index );
      if( index != nb_cls - 1 ) 
        bblas::noalias( bblas::subrange( temp, index, nb_cls ) )
           = bblas::subrange( *i_pis, index+1, nb_cls );
      *i_pis = temp; 
    }

    init_sums();
    construct_pairs();
      
    return index;
  }

  types::t_unsigned Regulated :: square_errors() const
  {
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( targets.size() != pis.size(),
              "Inconsistent number of returns and pis.\n" )

    types::t_real result(0);

    // loop over return values.
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.end();
    t_Weights :: const_iterator i_w = weights.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for(; i_target != i_target_end; ++i_target, ++i_w, ++i_pis)
    {
      __ASSERT( clusters.size() != i_pis->size(),
                "Inconsistent number of clusters and pis.\n" )

      types::t_real intermed( *i_target );
      // loop over alpha.
      t_StructPis :: const_iterator i_pi = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      t_Clusters :: const_iterator i_clusters = clusters.begin();
      for(; i_pi != i_pi_end; ++i_pi, ++i_clusters )
        intermed -= i_clusters->front().eci * (*i_pi);

      result += (*i_w) * intermed * intermed;
    } // end of loop over target values.

    return result / (types::t_real) targets.size(); 
  }

  void Regulated :: reassign( const t_Arg &_arg )
  {
    namespace bl = boost::lambda;

    // now copies to clusters.
    t_Clusters :: iterator i_clusters = clusters.begin();
    t_Clusters :: iterator i_clusters_end = clusters.end();
    t_Arg :: const_iterator i_arg = _arg.begin();
    for(; i_clusters != i_clusters_end; ++i_clusters, ++i_arg )
      std::for_each
      (
        i_clusters->begin(), i_clusters->end(),
        bl::bind( &Cluster::eci, bl::_1 ) = *i_arg
      );
  }

  types::t_real Regulated :: fit( t_Vector &_x ) 
  {
    namespace bl = boost::lambda;
    namespace bblas = boost::numeric::ublas;
    __ASSERT( nb_cls != clusters.size(),
              "Inconsistent number of clusters.\n" )
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of targets and weights.\n" )
    __ASSERT( pis.size() != targets.size(),
              "Inconsistent number of targets and pis.\n" )

    t_Matrix A = psums;
    cgs( psums, _x, esums );

    // computes square errors.
    types::t_real result(0), dummy;
    opt::concurrent_loop
    (
      weights.begin(), weights.end(), pis.begin(), targets.begin(),
      bl::var( result ) +=
      (
        bl::var( dummy ) = bl::_3 - bl::bind( &details::innerprod, bl::_2, _x ),
        bl::var( dummy ) * bl::var( dummy ) * bl::_1 
      )
    );
    return result / types::t_real( targets.size() );
  }

  void Regulated :: construct_pairs()
  {
    namespace bl = boost::lambda;
    if( fittingpairs.size() != targets.size() )
      fittingpairs.resize( targets.size() );
    types::t_unsigned k(0);
    opt::concurrent_loop
    (
       fittingpairs.begin(), fittingpairs.end(), k,
       bl::bind( &Regulated :: construct_pair, bl::var(this), bl::_1, bl::_2 )
    );
  }

  void Regulated :: construct_pair( t_FittingPair &_pair, types::t_unsigned &_k )
  {
    namespace bl = boost::lambda;
    __ASSERT( _k >= targets.size(),
              "Index out of range: " << _k << " >= " << targets.size() << ".\n" )

    // For easier identification
    t_FittingPair::first_type &_A = _pair.first;
    t_FittingPair::second_type &_b = _pair.second;

    // initializes matrix and array
    namespace bblas = boost::numeric::ublas;
    if( _A.size1() != nb_cls or _A.size2() != nb_cls )
      _A.resize( nb_cls, nb_cls );
    if( _b.size() != nb_cls ) _b.resize( nb_cls );

    // Values of removed structure.
    const types::t_real &target = targets[_k];
    const types::t_real &weight = weights[_k];
    const t_StructPis &kpis = pis[_k];
    
    // constructs b.
    _b = esums - weight * target * kpis;

    // loop over alpha.
    t_Matrix :: array_type :: const_iterator i_psum = psums.data().begin();
    t_StructPis :: const_iterator i_pi_begin = kpis.begin();
    t_StructPis :: const_iterator i_pi_end = kpis.end();
    t_StructPis :: const_iterator i_alphapi( i_pi_begin );
    t_Matrix :: array_type :: iterator i_A = _A.data().begin();
    for(; i_alphapi != i_pi_end; ++i_alphapi )
    {
      // loop over beta
      t_StructPis :: const_iterator i_betapi( i_pi_begin );
      for(; i_betapi != i_pi_end; ++i_A, ++i_betapi, ++i_psum )
      {
        __ASSERT( i_psum == psums.data().end(), "Iterator out of range.\n" )
        __ASSERT( i_A == _A.data().end(), "Iterator out of range.\n" )
        *i_A += (*i_psum) - weight * (*i_alphapi) * (*i_betapi );
      } // end of loop over beta
    } // end of loop over alpha.
  }

} // end of namespace CE

