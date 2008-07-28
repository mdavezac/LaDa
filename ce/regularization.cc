//
//  Version: $Id: functional_builder.impl.h 685 2008-07-22 17:20:39Z davezac $
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<algorithm>
#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>

#include <opt/random.h>
#include <opt/cgs.h>
#include <opt/gsl_nllsq.h>

#include "regularization.h"

namespace CE
{
  void Regulated :: init_sums()
  {
    namespace bl = boost::lambda;
    __ASSERT( pis.size() != targets.size(),
              "Inconsistent number of targets and pis.\n" )
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of targets and weights.\n" )

    // Resizes the clusters and fills with zeros.
    numsums.resize( nb_cls );
    denumsums.resize( nb_cls );
    std::fill( numsums.begin(), numsums.end(), 0e0 );
    t_DenumSums :: iterator i_denums = denumsums.begin();
    t_DenumSums :: iterator i_denums_end = denumsums.begin();
    for(; i_denums != i_denums_end; ++i_denums )
    {
      i_denums->resize( nb_cls );
      std::fill( i_denums->begin(), i_denums->end(), 0e0 );
    }

    // loop over targets.
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.end();
    t_Pis :: const_iterator i_pis = pis.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    for(; i_target != i_target_end; ++i_target, ++i_pis, ++i_w )
    {
      __ASSERT( i_pis->size() == nb_cls, "Inconsistent number of targets pis.\n" )
      // loop over alpha.
      t_NumSums :: iterator i_num = numsums.begin();
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      i_denums = denumsums.begin();
      for(; i_denums != i_denums_end; ++i_num, ++i_denums, ++i_alphapi )
      {
        *i_num += (*i_w) * (*i_target) * (*i_alphapi);

        // loop over betas.
        t_DenumSums :: value_type :: iterator i_denum = i_denums->begin();
        t_DenumSums :: value_type :: iterator i_denum_end = i_denums->end();
        t_StructPis :: const_iterator i_betapi = i_pis->begin();
        for(; i_denum != i_denum_end; ++i_denum, ++i_betapi )
          *i_denum += (*i_w) * (*i_alphapi) * (*i_betapi);

      } // end of loop over alphas.
    } // end of loop over betas.
  } 

  void Regulated :: operator()( const types::t_real * _arg, 
                                types::t_real *_return ) const
  {
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( targets.size() != pis.size(),
              "Inconsistent number of returns and pis.\n" )

    // loop over return values.
    types::t_real* i_ret = _return;
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.end();
    t_Targets :: const_iterator i_w = weights.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for(; i_target != i_target_end; ++i_ret, ++i_target, ++i_w, ++i_pis)
    {
      __ASSERT( i_pis->size() == numsums.size(),
                "Inconsistent number of pis and numsums.\n" )
      __ASSERT( i_pis->size() == denumsums.size(),
                "Inconsistent number of pis and denumsums.\n" )
      *i_ret = (*i_w) * (*i_target);

      // loop over alpha.
      const types::t_real * i_arg = _arg;
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      t_DenumSums :: const_iterator i_denums = denumsums.begin();
      for(; i_alphapi != i_pi_end; ++i_alphapi, ++i_denums )
      {
        // loop over beta.
        t_StructPis :: const_iterator i_betapi = i_pis->begin();
        t_DenumSums :: value_type :: const_iterator i_denum = i_denums->begin();
        t_NumSums :: const_iterator i_num = numsums.begin();
        for(; i_betapi != i_pi_end; ++i_betapi, ++i_num, ++i_denum )
        {
          types::t_real num, denum;
          num  = (*i_num) - (*i_w) * (*i_target) * (*i_betapi );
          num *= (*i_alphapi);
          denum = (*i_denum) - (*i_w) * (*i_alphapi) * (*i_betapi);
          if( i_alphapi == i_betapi ) denum += (*i_arg) * (*i_arg);
          
          *i_ret -= num / denum;
        } // end of loop over beta.
      }  // end of loop over alpha.
    } // end of loop over return values.
  }

  void Regulated :: jacobian( const types::t_real * _arg,
                              types::t_real *_jacobian ) const
  {
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( targets.size() != pis.size(),
              "Inconsistent number of returns and pis.\n" )

    // loop over return values.
    types::t_real *i_jac = _jacobian;
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.begin();
    t_Weights :: const_iterator i_w = weights.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for(; i_target != i_target_end; ++i_target, ++i_w, ++i_pis)
    {
      __ASSERT( i_pis->size() == numsums.size(),
                "Inconsistent number of pis and numsums.\n" )
      __ASSERT( i_pis->size() == denumsums.size(),
                "Inconsistent number of pis and denumsums.\n" )

      // loop over alpha.
      const types::t_real* i_arg = _arg;
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      t_DenumSums :: const_iterator i_denums = denumsums.begin();
      t_NumSums :: const_iterator i_num = numsums.begin();
      for( size_t i(0); i_alphapi != i_pi_end; 
           ++i_alphapi, ++i_denums, ++i_jac, ++i, ++i_num, ++i_arg )
      {
        __ASSERT( i_denums->size() <= i,
                  "Could not access fiagonal element of denumsums.\n" )
        types::t_real num, denum;
        num  = (*i_num) - (*i_w) * (*i_target) * (*i_alphapi );
        num *= (*i_alphapi);
        denum =   (*i_denums)[i] - (*i_w) * (*i_alphapi) * (*i_alphapi)
                + (*i_arg) * (*i_arg);
        
        *i_jac = num / ( denum * denum );
      }  // end of loop over alpha.
    } // end of loop over targets.
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

    nb_cls = clusters.size();
  } 
 
  types::t_unsigned Regulated :: reduce()
  {
     typedef std::vector< std::vector<Cluster> > t_Clusters;
     t_Clusters :: iterator i_cls = clusters.begin(); 
     t_Clusters :: iterator i_cls_end = clusters.end(); 
     t_Clusters :: iterator i_found( i_cls ); ++i_cls;
     types::t_unsigned index(0);
     types::t_real least( std::abs(i_cls->front().eci ) );
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

       Regulated :: t_StructPis :: iterator i_find = i_pis->begin();
       for( types::t_unsigned i(index); i > 0; --i, ++i_find );
       i_pis->erase( i_find );
     }

     __ASSERT( numsums.size() != denumsums.size(),
               "Inconsistent number of numsums and denumsums.\n" )
     __ASSERT( numsums.size() <= index, "Index out of range.\n" )
     t_NumSums :: iterator i_num = numsums.begin();
     t_DenumSums :: iterator i_denums = denumsums.begin();
     for(size_t i(index); i > 0; ++i_num, ++i_denums, --i );
     numsums.erase( i_num );
     denumsums.erase( i_denums );
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
    t_Targets :: const_iterator i_target_end = targets.begin();
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
        intermed -= i_clusters->size() * i_clusters->front().eci * (*i_pi);

      result += (*i_w) * intermed * intermed;
    } // end of loop over target values.

    return result / (types::t_real) targets.size(); 
  }

  void Regulated :: reassign( const t_Arg &_arg )
  {
    namespace bl = boost::lambda;
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of returns and weights.\n" )
    __ASSERT( targets.size() != pis.size(),
              "Inconsistent number of returns and pis.\n" )
    __ASSERT( _arg.size() != clusters.size(),
              "Inconsistent number of returns and weights.\n" )

    // Since clusters can/are lists, we need to jump through a couple of hoops.
    std::vector< types::t_real > ecis( nb_cls, 0e0 );

    // loop over target values.
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.end();
    t_Targets :: const_iterator i_w = weights.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for(; i_target != i_target_end; ++i_target, ++i_w, ++i_pis)
    {
      __ASSERT( i_pis->size() == numsums.size(),
                "Inconsistent number of pis and numsums.\n" )
      __ASSERT( i_pis->size() == denumsums.size(),
                "Inconsistent number of pis and denumsums.\n" )
      // loop over alpha.
      t_Arg :: const_iterator i_arg = _arg.begin();
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      t_DenumSums :: const_iterator i_denums = denumsums.begin();
      std::vector< types::t_real > :: iterator i_eci = ecis.begin();
      for(; i_alphapi != i_pi_end; ++i_alphapi, ++i_denums, ++i_eci )
      {
        // loop over beta.
        t_StructPis :: const_iterator i_betapi = i_pis->begin();
        t_DenumSums :: value_type :: const_iterator i_denum = i_denums->begin();
        t_NumSums :: const_iterator i_num = numsums.begin();
        for(; i_betapi != i_pi_end; ++i_betapi, ++i_num, ++i_denum )
        {
          types::t_real num, denum;
          num  = (*i_num) - (*i_w) * (*i_target) * (*i_betapi );
          denum = (*i_denum) - (*i_w) * (*i_alphapi) * (*i_betapi);
          if( i_alphapi == i_betapi ) denum += (*i_arg) * (*i_arg);
          
          *i_eci += num / denum;
        } // end of loop over beta.
      }  // end of loop over alpha.
    } // end of loop over return values.

    // now copies to clusters.
    t_Clusters :: iterator i_clusters = clusters.begin();
    t_Clusters :: iterator i_clusters_end = clusters.end();
    std::vector< types::t_real > :: const_iterator i_eci = ecis.begin();
    for(; i_clusters != i_clusters_end; ++i_clusters, ++i_eci )
      std::for_each
      (
        i_clusters->begin(), i_clusters->end(),
          bl::bind( &Cluster::eci, bl::_1 )
        = bl::constant( *i_eci / (types::t_real) i_clusters->size() )
      );
  }


  types::t_real Regulated :: fit( types::t_real _tol,
                                  types::t_unsigned _imax,
                                  bool _verbose ) 
  {
    namespace bl = boost::lambda;
    __ASSERT( nb_cls != clusters.size(),
              "Inconsistent number of clusters.\n" )
    __ASSERT( targets.size() != weights.size(),
              "Inconsistent number of targets and weights.\n" )
    __ASSERT( pis.size() != targets.size(),
              "Inconsistent number of targets and pis.\n" )

    // initializes matrix and array
    namespace bblas = boost::numeric::ublas;
    bblas :: matrix<types::t_real> A( nb_cls, nb_cls );
    bblas :: vector<types::t_real> b( nb_cls );
    bblas :: vector<types::t_real> x( nb_cls );
    std::fill( A.data().begin(), A.data().end(), 0e0 );
    std::fill( b.data().begin(), b.data().end(), 0e0 );

    std::vector< types::t_real > multip( nb_cls );
    std::transform
    (
      clusters.begin(), clusters.end(), multip.begin(),
      bl::bind<types::t_real>( &t_Clusters::value_type::size, bl::_1 ) 
    );

    // loop over targets.
    t_Targets :: const_iterator i_target = targets.begin();
    t_Targets :: const_iterator i_target_end = targets.end();
    t_Weights :: const_iterator i_w = weights.begin();
    t_Pis :: const_iterator i_pis = pis.begin();
    for(; i_target != i_target_end; ++i_target, ++i_w, ++i_pis )
    {
      __ASSERT( i_pis->size() != nb_cls,
                "Inconsistent number of pis and clusters.\n" )
      // loop over alpha
      bblas::vector<types::t_real> :: iterator i_b = b.begin();
      std::vector< types::t_real > :: const_iterator i_alphamul = multip.begin();
      t_StructPis :: const_iterator i_alphapi = i_pis->begin();
      t_StructPis :: const_iterator i_pi_end = i_pis->end();
      for( size_t alpha(0); i_alphapi != i_pi_end; 
           ++i_b, ++i_alphapi, ++i_alphamul, ++alpha )
      {
        // adds to target vector.
        types::t_real row( (*i_w) * (*i_target) * (*i_alphamul ) * (*i_alphapi) );
        *i_b += row;

        // loop over beta
        std::vector< types::t_real > :: const_iterator i_betamul = multip.begin();
        t_StructPis :: const_iterator i_betapi = i_pis->begin();
        for( size_t beta(0); i_betapi != i_pi_end; ++i_betapi, ++i_betamul, ++beta )
          A( alpha, beta ) += row * (*i_betamul) * (*i_betapi );
      } // end of loop over alpha
    } // end of loop over structures.

    // Now initialize and run conjugate gradient procedure.
    Fitting::Cgs cgs;
    cgs.tolerance = _tol;
    cgs.itermax = _imax;
    cgs.verbose = _verbose;
    cgs( A, x, b );

    // reassigns ecis.
    bblas::vector<types::t_real>::array_type :: const_iterator
      i_x = x.data().begin();
    bblas::vector<types::t_real>::array_type :: const_iterator
      i_x_end = x.data().end();
    t_Clusters :: iterator i_clusters = clusters.begin();
    for(; i_x != i_x_end; ++i_x )
      std::for_each
      (
        i_clusters->begin(), i_clusters->end(),
          bl::bind( &t_Clusters::value_type::value_type::eci, bl::_1 ) 
        = bl::constant( *i_x ) 
      );

    // returns error.
    return square_errors();
  }

  void variational_figure_plot( Regulated &_reg, 
                                types::t_real _tolerance,
                                types::t_int _verbosity )
  {
    namespace bl = boost::lambda;

    std::cout << "Starting Variational figure plot.\n";

    Regulated::t_Clusters save_clusters( _reg.clusters );

    Fitting :: NonLinearGsl solver;
    solver.tolerance = _tolerance;
    solver.verbose = _verbosity >= 2;
    
    types::t_unsigned nb_cls( _reg.clusters.size() );
    Fitting :: NonLinearGsl:: t_Vector solution( nb_cls );
    const Fitting :: NonLinearGsl :: t_Vector zero_vec( solution.size(), 0e0);
    Fitting :: NonLinearGsl :: t_Vector cv_vec( _reg.dim() );
    const types::t_real range(1);
    std::for_each
    ( 
      solution.begin(), solution.end(),
      bl::_1 =   bl::bind( &opt::random::rng ) * bl::constant(range)
               + bl::constant( 1e0 - 0.5*range )
    );

    while( _reg.clusters.size() > 0 )
    {
      // Fitting Error
      types::t_real fit = _reg.fit( _tolerance, 40, _verbosity >= 3);

      // CV with zero weights
      _reg( &zero_vec[0], &cv_vec[0]);
      types::t_real cvwz(0); 
      std::for_each
      ( 
        cv_vec.begin(), cv_vec.end(),
        bl::var(cvwz) += bl::_1 * bl::_1
      );

      // CV with optimized weights
      types::t_real cvw = solver( _reg, solution );
      
      // reduces figures by one.
      _reg.reassign( solution );
      types::t_unsigned index = _reg.reduce();

      __ASSERT( index > nb_cls, "index out-of-range.\n" )
      Regulated :: t_Clusters :: const_iterator i_found( save_clusters.begin() );
      for( size_t i(index); i > 0; --i, ++i_found );

      std::cout << " Number of clusters: " << index << "\n"
                << "  CV with weights: " << cvw << "\n"
                << "  CV with weights=0: " << cvwz << "\n"
                << "  Fitting squared error: " << fit << "\n"
                << "  Dropping cluster " << index << "\n"
                << i_found->front() << "\n\n";
    }
  }

} // end of namespace CE

