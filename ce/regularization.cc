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
    __DEBUGTRYBEGIN
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
      for(size_t i(0); i_alphapi != i_pi_end; ++i_alphapi, ++i_esum, ++i)
      {
        __ASSERT( i_esum == esums.end(), "Iterator out of range.\n" )
        esums(i) += (*i_w) * (*i_target) * (*i_alphapi);

        // loop over betas.
        t_StructPis :: const_iterator i_betapi( i_pi_begin );
        for(size_t j(0); i_betapi != i_pi_end; ++i_betapi, ++i_psum, ++j )
        {
          __ASSERT( i_psum == psums.data().end(), "Iterator out of range.\n" )
          psums(i,j) += (*i_w) * (*i_alphapi) * (*i_betapi);
        }
      } // end of loop over alphas.
    } // end of loop over betas.
    __DEBUGTRYEND(, "Error in Regulated::init_sums.\n" )
  } 

  Regulated::t_Return Regulated :: operator()( const types::t_real * _arg ) const
  {
    __DEBUGTRYBEGIN
    namespace bblas = boost::numeric::ublas;
    namespace bl = boost::lambda;
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( targets.size() != fittingpairs.size(),
              "Incorrect number of training pairs.\n" )

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

  void Regulated :: init( const t_Structures &_structures )
  {
    __DEBUGTRYBEGIN
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
    __DEBUGTRYEND(, "Error in Regulated::init().\n" )
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

    t_Pis :: iterator i_pis = pis.begin();
    t_Pis :: iterator i_pis_end = pis.end();
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

    init_sums();
    construct_pairs();
      
    return result;
    __DEBUGTRYEND(, "Error in Regulated::reduce().\n" )
  }

  types::t_unsigned Regulated :: square_errors() const
  {
    __DEBUGTRYBEGIN
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
    __DEBUGTRYEND(, "Error in Regulated::square_errors().\n" )
  }

  void Regulated :: reassign( const t_Arg &_arg )
  {
    __DEBUGTRYBEGIN
    namespace bl = boost::lambda;

    // now copies to clusters.
    t_Clusters :: iterator i_clusters = clusters.begin();
    t_Clusters :: iterator i_clusters_end = clusters.end();
    t_Arg :: const_iterator i_arg = _arg.begin();
    for(; i_clusters != i_clusters_end; ++i_clusters, ++i_arg )
      std::for_each
      (
        i_clusters->begin(), i_clusters->end(),
        bl::bind( &Cluster::eci, bl::_1 ) = bl::constant(*i_arg)
      );
    __DEBUGTRYEND(, "Error in Regulated::reassign().\n" )
  }

  types::t_real Regulated :: fit( t_Vector &_x, const t_Vector &_weights ) const
  {
    __DEBUGTRYBEGIN
    namespace bl = boost::lambda;
    namespace bblas = boost::numeric::ublas;
    __ASSERT( nb_cls != clusters.size(),
              "Inconsistent number of clusters.\n" )
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of targets and weights.\n" )
    __ASSERT( pis.size() != targets.size(),
              "Inconsistent number of targets and pis.\n" )

    t_Matrix A = psums;
    Regulated::t_Vector :: const_iterator i_arg = _weights.begin();
    Regulated::t_Vector :: const_iterator i_arg_end = _weights.end();
    for( size_t i(0); i_arg != i_arg_end; ++i_arg, ++i )
      A(i,i) += (*i_arg) * (*i_arg);
    cgs( A, _x, esums );

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
    __DEBUGTRYEND(, "Error in Regulated::fit().\n" )
  }

  void Regulated :: construct_pairs()
  {
    __DEBUGTRYBEGIN
    namespace bl = boost::lambda;
    if( fittingpairs.size() != targets.size() )
      fittingpairs.resize( targets.size() );
    ;
    t_FittingPairs :: iterator i_fpair = fittingpairs.begin();
    t_FittingPairs :: iterator i_fpair_end = fittingpairs.end();
    for(types::t_unsigned k(0); i_fpair != i_fpair_end; ++i_fpair, ++k )
      construct_pair( *i_fpair, k );
    __DEBUGTRYEND(, "Error in Regulated::construct_pairs().\n" )
  }

  void Regulated :: construct_pair( t_FittingPair &_pair, types::t_unsigned &_k )
  {
    __DEBUGTRYBEGIN
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
        *i_A = (*i_psum) - weight * (*i_alphapi) * (*i_betapi );
      } // end of loop over beta
    } // end of loop over alpha.
    __DEBUGTRYEND(, "Error in Regulated::construct_pair().\n" )
  }

  void Fit :: leave_one_out( const Regulated &_reg )
  {
    __DEBUGTRYBEGIN
    namespace bblas = boost::numeric::ublas;
    namespace bl = boost::lambda;

    opt::ErrorTuple training(0,0,0), prediction(0,0,0);

    std::cout << "\nLeave-one-out procedure\n";

    // computes mean and variance.
    types::t_real norm( 0 );
    opt::NErrorTuple nerror( opt::mean_n_var( _reg.targets, _reg.weights ) );
    std::cout << "  Weighted mean of data: " << nerror.mean << "\n" 
              << "  Weighted variance of data: " << nerror.variance << "\n";


    Regulated :: t_Vector ecis(_reg.nb_cls );
    Regulated :: t_Vector weights(_reg.nb_cls );
    std::fill( weights.begin(), weights.end(), 0e0 );
    std::fill( ecis.begin(), ecis.end(), 0e0 );
    if( do_pairreg ) pair_reg( _reg, weights );

    // loop over structures
    Regulated :: t_Weights :: const_iterator i_w = _reg.weights.begin();
    for( types::t_int n = 0; n < structures.size(); ++n, ++i_w )
    {
      fit_but_one( _reg, ecis, weights, n );
      if( verbose ) std::cout <<  "Training errors:\n" ;
      opt::ErrorTuple train( check_all( _reg, ecis, n ) );
      training += train;
      if( verbose ) std::cout << "    " << ( nerror = train ) << "\n"; 
      if( verbose ) std::cout <<  "Prediction errors:\n" ;
      opt::ErrorTuple pred( check_one( _reg, ecis, n ), *i_w );
      if( verbose ) std::cout << "    " << ( nerror = pred ) << "\n"; 
      prediction += pred;
      norm += *i_w;
    }
    prediction.get<0>() /= norm;
    prediction.get<1>() /= norm;
    training.get<0>() /= norm;
    training.get<1>() /= norm;
    std::cout << "\nTraining Errors:\n" << "   " << ( nerror = training ) << "\n";
    std::cout << "Prediction Errors:\n" << "   " << ( nerror = prediction ) << "\n";
    __DEBUGTRYEND(, "Error in Fit::leave_one_out().\n" )
  }
  void Fit :: fit( const Regulated &_reg )
  {
    __DEBUGTRYBEGIN
    namespace bblas = boost::numeric::ublas;
    namespace bl = boost::lambda;

    std::cout << "\nFitting procedure\n";
    opt::NErrorTuple nerror( opt::mean_n_var( _reg.targets, _reg.weights ) );
    std::cout << "  Weighted mean of data: " << nerror.mean << "\n" 
              << "  Weighted variance of data: " << nerror.variance << "\n";

    Regulated :: t_Vector ecis(_reg.nb_cls );
    Regulated :: t_Vector weights(_reg.nb_cls );
    std::fill( weights.begin(), weights.end(), 0e0 );
    std::fill( ecis.begin(), ecis.end(), 0e0 );
    if( do_pairreg ) pair_reg( _reg, weights );

    // fits.
    _reg.fit( ecis, weights );
    std::cout << "\nTraining Errors:\n";
    opt::ErrorTuple training( check_all( _reg, ecis ) );
    std::cout << ( nerror = training ) << "\n";
    __DEBUGTRYEND(, "Error in Fit::fit().\n" )
  }

  void Fit :: pair_reg( const Regulated &_reg, Regulated :: t_Vector &_weights )
  {
    __DEBUGTRYBEGIN
    Regulated :: t_Clusters :: const_iterator i_clusters = _reg.clusters.begin();
    Regulated :: t_Clusters :: const_iterator i_clusters_end = _reg.clusters.end();
    Regulated :: t_Vector :: iterator i_weight = _weights.begin();
     
    // Now loop over pairs.
    types::t_real normalization(0);
    for(; i_clusters != i_clusters_end; ++i_clusters, ++i_weight )
    {
      __ASSERT( i_clusters->size() == 0, "No clusters in class.\n" )
      if( i_clusters->front().size() != 2 ) continue;

      types::t_real D( i_clusters->size() );
      types::t_real R( atat::norm2(   i_clusters->front()[0] 
                                    - i_clusters->front()[1] ) );

      *i_weight = std::pow( R, lambda ) / D;
      normalization += laksreg ?
                         ( std::pow( R, lambda ) * D ):
                         ( std::sqrt( std::pow( R, lambda  ) / D ) );
      *i_weight = std::sqrt( *i_weight );
    }
    normalization = laksreg ?
                      std::sqrt( tcoef / normalization ):
                      std::sqrt( tcoef ) / normalization;
    i_clusters = _reg.clusters.begin();
    i_weight = _weights.begin();
    for(; i_clusters != i_clusters_end; ++i_clusters, ++i_weight )
    {
      if( i_clusters->front().size() != 2 ) continue;
      *i_weight *= normalization;
    }
    __DEBUGTRYEND(, "Error in Fit::pair_reg().\n" )
  }

  void Fit :: fit_but_one( const Regulated &_reg,
                           Regulated :: t_Vector &_x,
                           const Regulated :: t_Vector &_weights,
                           const types::t_unsigned _n ) const 
  {
    __DEBUGTRYBEGIN
    // loop over target values.
    const Regulated :: t_FittingPairs :: value_type pair = _reg.fittingpairs[_n];

    // construct matrix with weights.
    Regulated :: t_Matrix A( pair.first );
    Regulated::t_Vector :: const_iterator i_arg = _weights.begin();
    Regulated::t_Vector :: const_iterator i_arg_end = _weights.end();
    for( size_t i(0); i_arg != i_arg_end; ++i_arg, ++i )
      A(i,i) += (*i_arg) * (*i_arg);

    // fits.
    _reg.cgs( A, _x, pair.second );
    __DEBUGTRYEND(, "Error in Fit::fit_but_one().\n" )
  }

  types::t_real Fit :: check_one( const Regulated &_reg, 
                                  const Regulated :: t_Vector &_ecis,
                                  types::t_unsigned _n )
  {
    __DEBUGTRYBEGIN
    namespace fs = boost::filesystem;
    namespace bblas = boost::numeric::ublas;
    const Regulated :: t_Targets :: value_type &target = _reg.targets[_n];
    const Regulated :: t_Weights :: value_type &weight = _reg.weights[_n];
    const Regulated :: t_Pis :: value_type pis = _reg.pis[_n];
    const std::string name = fs::path( structures[_n].name ).leaf() ;
 
    types::t_real predic( bblas::inner_prod( _ecis, pis ) );
    types::t_real error( std::abs( target - predic ) );
    if( verbose )
      std::cout << "  structure: " << std::setw(30) << name << "  "
                << "Target: " << std::fixed << std::setw(8) 
                << std::setprecision(2) << target << " "
                << "Separable: " << std::fixed << std::setw(8)
                << std::setprecision(2) << predic << "   "
                << "|Target-Separable| * weight: "
                << std::fixed << std::setw(10) << std::setprecision(3) 
                << weight * error
                << "\n";
    return error;
    __DEBUGTRYEND(, "Error in Fit::check_one().\n" )
  }

  opt::ErrorTuple Fit :: check_all( const Regulated &_reg, 
                                    const Regulated :: t_Vector &_ecis,
                                    types::t_int _n  )
  {
    __DEBUGTRYBEGIN
    opt::ErrorTuple result( 0,0,0 );
    types::t_real norm(0);
    Regulated :: t_Weights :: const_iterator i_weight = _reg.weights.begin();
    for(types::t_int n(0); n < structures.size(); ++n, ++i_weight )
    {
      if( _n == n ) continue;
      types::t_real error( check_one( _reg, _ecis, n ) );
      result += opt::ErrorTuple( error, (*i_weight) );
      norm += (*i_weight);
    }
    types::t_real N( structures.size() );
    result.get<0>() /= norm;
    result.get<1>() /= norm;
    return result;
    __DEBUGTRYEND(, "Error in Fit::check_all().\n" )
  }

} // end of namespace CE

