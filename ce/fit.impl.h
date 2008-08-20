//
//  Version: $Id$
//

#include <boost/lambda/bind.hpp>
#include <opt/random.h>

namespace CE
{

  namespace details
  {
    template< class T_SOLVER, class T_CLASS >
      opt::ErrorTuple operator_( const T_CLASS &_class,
                                 BaseFit::t_Vector &_x, 
                                 const T_SOLVER &_solver ) 
      {
        __DEBUGTRYBEGIN
        BaseFit::t_Matrix A( _class.nb_cls, _class.nb_cls );
        BaseFit::t_Vector b( _class.nb_cls );
        std::fill( b.begin(), b.end(), 0e0 );
        std::fill( A.data().begin(), A.data().end(), 0e0 );
        _class.construct( A, b );
        _solver( A, _x, b );
  
        // computes square errors.
        opt::ErrorTuple result = _class.check_training( _x, _class.verbose );
        return result;
        __DEBUGTRYEND(, "Error in CE::details::operator_().\n" )
      }
    template< class T_CLASS, class T_MATRIX, class T_VECTOR >
      void construct_( const T_CLASS &_class, T_MATRIX &_A, T_VECTOR &_b ) 
      {
        __DEBUGTRYBEGIN
        __ASSERT( _A.size1() != _class.nb_cls,
                     "Inconsistent sizes: " << _A.size1()
                  << " != " << _class.nb_cls << "\n" )
        __ASSERT( _A.size2() != _class.nb_cls, 
                     "Inconsistent sizes: " << _A.size2()
                  << " != " << _class.nb_cls << "\n" )
        __ASSERT( _b.size() != _class.nb_cls,
                     "Inconsistent sizes: " << _b.size()
                  << " != " << _class.nb_cls << "\n" )
        _class.create_A_n_b( _A, _b );
        _class.other_A_n_b( _A, _b );
        __DEBUGTRYEND(, "Error in CE::details::construct_().\n" )
      }
  }

  template<class T_POLICY> template< class T_MATRIX, class T_VECTOR>
  void Fit<T_POLICY> :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b ) const
  {
    __DEBUGTRYBEGIN
    namespace bl = boost::lambda;
    namespace bblas = boost::numeric::ublas;
    __ASSERT( pis.size() != structures().size(),
              "Inconsistent number of structures and pis.\n" )
    __ASSERT( structures().size() != weights().size(),
              "Inconsistent number of structures and weights.\n" )

    // loop over targets.
    t_Structures :: const_iterator i_target = structures().begin();
    t_Structures :: const_iterator i_target_end = structures().end();
    t_Pis :: const_iterator i_pis = pis.begin();
    t_Weights :: const_iterator i_w = weights().begin();
    for( types::t_unsigned i(0); 
         i_target != i_target_end; 
         ++i_target, ++i_pis, ++i_w, ++i )
    {
      if( t_Policy :: found( i ) ) continue; 
      _A += (*i_w) * bblas::outer_prod( *i_pis, *i_pis );
      _b += (*i_w) * i_target->energy * (*i_pis);
    } // end of loop over betas.
    __DEBUGTRYEND(, "Error in Fit::create_A_n_b.\n" )
  } 


  template< class T_BASE > template< class T_SOLVER >
   opt::ErrorTuple RegulatedFit<T_BASE>
     :: operator()( BaseFit::t_Vector &_x,
                    const types::t_real *_weights,
                    T_SOLVER &_solver ) const
     {
       regweights = _weights;
       opt::ErrorTuple result = details::operator_( *this, _x, _solver );
       regweights = NULL;
       return result;
     }
  template< class T_BASE > template< class T_MATRIX, class T_VECTOR >
    void RegulatedFit<T_BASE> :: other_A_n_b( T_MATRIX &_A, T_VECTOR &_b ) const
    {
      t_Base :: other_A_n_b( _A, _b );
      if( not regweights  ) return;
      __DEBUGTRYBEGIN
      const types::t_real *i_w = regweights;
      const types::t_real *i_w_end = regweights + nb_cls;
      for(size_t i(0); i_w != i_w_end; ++i_w, ++i )
        _A( i, i ) += (*i_w) * (*i_w);
      __DEBUGTRYEND(, "Error in PairReg<T_POLICY>::other_A_n_b().\n" )
    }

  namespace FittingPolicy
  {
    template< class T_BASE >
      bool Excluded<T_BASE> :: found( types::t_unsigned _i ) const 
      {
        if( not excluded.size() ) return false;
        return std::find( excluded.begin(), excluded.end(), _i ) != excluded.end();
      }
    template< class T_BASE >
      opt::ErrorTuple Excluded<T_BASE> :: check( const BaseFit::t_Vector &_ecis,
                                                 bool _training, bool _verbose ) const
      {
        __DEBUGTRYBEGIN
        opt::ErrorTuple result;
        for(types::t_int n(0); n < structures().size(); ++n)
          if( found(n) xor _training ) 
            result += t_Base::check_one( _ecis, n, _verbose );
        return result;
        __DEBUGTRYEND(, "Error in Fit::check_all().\n" )
      }

    template< class T_BASE >
    opt::ErrorTuple Nothing<T_BASE> :: check_training( const BaseFit::t_Vector &_ecis,
                                                       bool _verbose )
    {
      __DEBUGTRYBEGIN
      opt::ErrorTuple result;
      for(types::t_int n(0); n < structures().size(); ++n)
        result += t_Base::check_one( _ecis, n, _verbose );
      return result;
      __DEBUGTRYEND(, "Error in Fit::check_all().\n" )
    }

    template< class T_BASE >
      void PairReg<T_BASE> :: init( const BaseFit::t_Clusters &_clusters )
      {
        t_Base::init( _clusters );
        __DEBUGTRYBEGIN
        BaseFit::t_Clusters::const_iterator i_clusters = _clusters.begin();
        BaseFit::t_Clusters::const_iterator i_clusters_end = _clusters.end();
         
        // Now loop over pairs.
        pairweights.clear();
        types::t_real normalization(0);
        for(size_t i(0); i_clusters != i_clusters_end; ++i_clusters, ++i )
        {
          __ASSERT( i_clusters->size() == 0, "No clusters in class.\n" )
          if( i_clusters->front().size() != 2 ) continue;
  
          types::t_real D( i_clusters->size() );
          types::t_real R( atat::norm2(   i_clusters->front()[0] 
                                        - i_clusters->front()[1] ) );
  
          types::t_real weight = std::pow( R, alpha ) / D;
          normalization += laksreg ?
                             ( std::pow( R, alpha ) * D ):
                             ( std::sqrt( std::pow( R, alpha  ) / D ) );
          weight = std::sqrt( weight );
          pairweights.push_back( std::make_pair( i, weight ) );
        }
        normalization = laksreg ?
                          std::sqrt( tcoef / normalization ):
                          std::sqrt( tcoef ) / normalization;
        t_PairWeights :: iterator i_w = pairweights.begin();
        t_PairWeights :: iterator i_w_end = pairweights.end();
        for(; i_w != i_w_end; ++i_w ) i_w->second *= normalization;
        __DEBUGTRYEND(, "Error in PairReg<T_POLICY>::init().\n" )
      }

    template< class T_BASE > template< class T_MATRIX, class T_VECTOR >
      void PairReg<T_BASE> ::  other_A_n_b( T_MATRIX &_A, T_VECTOR &_b ) const 
      {
        t_Base :: other_A_n_b( _A, _b );
        __DEBUGTRYBEGIN
        if( not do_pairreg ) return;
        t_PairWeights :: const_iterator i_w = pairweights.begin();
        t_PairWeights :: const_iterator i_w_end = pairweights.end();
        for(; i_w != i_w_end; ++i_w )
        {
          __ASSERT( i_w->first >= _A.size1(),
                    "Inconsistent sizes: " << i_w->first << " >= " << _A.size1() << "\n" )
          __ASSERT( i_w->first >= _A.size2(),
                    "Inconsistent sizes: " << i_w->first << " >= " << _A.size2() << "\n" )
          _A( i_w->first, i_w->first ) += i_w->second * i_w->second;
        }
        __DEBUGTRYEND(, "Error in PairReg<T_POLICY>::other_A_n_b().\n" )
      }
  } // end of namespace FittingPolicy

  template< class T_FIT, class T_SOLVER >
    std::pair< opt::ErrorTuple, opt::ErrorTuple >
      leave_one_out( const T_FIT &_fit, const T_SOLVER &_solver,
                     BaseFit::t_Vector &_x, bool _verbose )
      {
        opt::ErrorTuple training, prediction;
        types::t_unsigned N = _fit.structures().size();
     
        _fit.excluded.resize(1);
        for(size_t i(0); i < N; ++i )
        {
          _fit.excluded[0] = i;
          opt::ErrorTuple intermediate;
          if( _verbose ) std::cout << "  " << _fit.excluded[0]
                                   << ". Training Errors: ";
          intermediate = _fit( _x, _solver );
          if( _verbose ) std::cout << intermediate << "\n";
          training += intermediate;
     
          if( _verbose ) std::cout << "  " << _fit.excluded[0]
                                   << ". Prediction Errors: ";
          intermediate = _fit.check_predictions( _x, _fit.verbose );
          if( _verbose ) std::cout << intermediate << "\n";
          prediction += intermediate;
        }
     
        return std::make_pair( training, prediction );
      }
  template< class T_FIT, class T_SOLVER >
    opt::t_ErrorPair leave_many_out( Fitting::LeaveManyOut &_lmo,
                                     const T_FIT &_fit, const T_SOLVER &_solver )
    {
      namespace bl = boost::lambda;
      __TRYBEGIN
      BaseFit::t_Vector x( _fit.dof() );
      opt::t_ErrorPair errors;
      if( not _lmo.do_perform ) return errors;
      if( not _lmo.sets.size() ) _lmo.create_sets( _fit.structures().size() );
 
      typedef std::vector< std::vector< types::t_unsigned > >
                                    :: const_iterator const_iterator;
      const_iterator i_set = _lmo.sets.begin();
      const_iterator i_set_end = _lmo.sets.end();
      for(size_t n(0); i_set != i_set_end; ++i_set, ++n )
      {
        const types::t_real range(10);
        std::for_each
        ( 
          x.begin(), x.end(),
          bl::_1 =   ( bl::bind( &opt::random::rng ) - bl::constant(5e-1) )
                   * bl::constant(range)
        );
        { // fitting.
          opt::ErrorTuple intermediate;
          _fit.excluded = *i_set;
          if( _lmo.verbosity >= 1 ) std::cout << " " << n
                                              << ". Training Errors: ";
          if( _lmo.verbosity >= 2 ) std::cout << "\n";
          intermediate = _fit( x, _solver );
          if( _lmo.verbosity >=1 ) std::cout << intermediate << "\n";
          errors.first += intermediate;
        } // end of fitting.
        { // Prediction.
          opt::ErrorTuple intermediate;
          if( _lmo.verbosity >= 1 ) std::cout << " " << n
                                          << ". Training Errors: ";
          if( _lmo.verbosity >= 2 ) std::cout << "\n";
          intermediate = _fit.check_predictions( x, _fit.verbose );
          if( _lmo.verbosity >=1 ) std::cout << intermediate << "\n";
          errors.second += intermediate;
        }
      }
      return errors;
      __TRYEND(, "Error while performing leave-many-out.\n")
    }
  template< class T_FIT, class T_SOLVER >
    std::pair< opt::ErrorTuple, opt::ErrorTuple >
      leave_one_out( const T_FIT &_fit,
                     const T_SOLVER &_solver,
                     BaseFit::t_Vector &_x, 
                     const types::t_real *_w, 
                     bool _verbose )
      {
        opt::ErrorTuple training, prediction;
        types::t_unsigned N = _fit.structures().size();
     
        _fit.excluded.resize(1);
        for(size_t i(0); i < N; ++i )
        {
          _fit.excluded[0] = i;
          opt::ErrorTuple intermediate;
          if( _verbose ) std::cout << "  " << _fit.excluded[0]
                                   << ". Training Errors: ";
          intermediate = _fit( _x, _w, _solver );
          if( _verbose ) std::cout << intermediate << "\n";
          training += intermediate;
     
          if( _verbose ) std::cout << "  " << _fit.excluded[0]
                                   << ". Prediction Errors: ";
          intermediate = _fit.check_predictions( _x, _fit.verbose );
          if( _verbose ) std::cout << intermediate << "\n";
          prediction += intermediate;
        }
     
        return std::make_pair( training, prediction );
      }

} // end of namespace CE
