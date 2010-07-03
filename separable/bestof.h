#ifndef _SEPARABLES_BESTOF_H_
#define _SEPARABLES_BESTOF_H_

#include "LaDaConfig.h"

#include <utility> 

#include <opt/types.h>
#include <math/fuzzy.h>
#include <opt/debug.h>


namespace LaDa
{
  namespace Separables
  {
    //! \brief Aggregates a number of results, and chooses the best.
    //! \tparam T_SOLVER A linear solver type.
    template< class T_SOLVER >
    class BestOf : public T_SOLVER
    {
      public:
        //! Type of the Matrix for collapsed 1d problem.
        typedef T_SOLVER t_Solver;

      public:
        //! Maximum number of iterations.
        types::t_unsigned n;
        //! Verbosity
        bool verbose;
        //! Wether to perform reruns as preruns.
        bool prerun;

      public:
        //! Constructor.
        BestOf   ( types::t_unsigned _n = 1 ) 
               : t_Solver(), n(_n), verbose(false), prerun(false) {}
        //! Copy Constructor.
        BestOf   ( BestOf &_c )
               : t_Solver( _c ), n(_c.n), verbose(_c.verbose), prerun(_c.prerun) {}
        //! Destructor
        ~BestOf() {}
        //! The functor itself.
        template< class T_COEFS, class T_COLLAPSE >
          types::t_real operator()( T_COEFS &_coefs, T_COLLAPSE *_collapse );
    };
    template< class T_SOLVER >
      template< class T_COEFS, class T_COLLAPSE > types::t_real
        BestOf<T_SOLVER> :: operator()( T_COEFS &_coefs, T_COLLAPSE *_collapse )
        {
          __DOASSERT( n==0, "0 runs required on input. We are done here.\n" )
          if( n == 1 ) return t_Solver :: operator()( _coefs, _collapse );

          bool docopy( false );
          typedef typename std::pair< types::t_real, T_COEFS > t_Best;
          t_Best best( t_Solver :: operator()( _coefs, _collapse ),
                       _coefs );
          typename T_COLLAPSE::t_Function::t_Coefs
            norms( _collapse->function.coefs );
          size_t fieldsize = (size_t) std::ceil( std::log( (types::t_real) n)/std::log(10e0) );
          ++fieldsize;
          if( verbose ) std::cout << "Run " << std::setw(fieldsize)
                                  << 1 << ": " << std::setw(15)
                                  << std::fixed
                                  << std::setprecision(3) << best.first << "\n";
          for( types::t_unsigned i(2); i <= n; ++i )
          {
            docopy = true;
            _collapse->create_coefs( _coefs );
            types::t_real result( t_Solver ::operator()( _coefs, _collapse ) );

            if( verbose ) std::cout << "Run " << std::setw(fieldsize)
                                    << i << ": " << std::setw(15)
                                    << std::fixed
                                    << std::setprecision(3) << result << "\n";

            if( math::geq( result, best.first ) ) continue;
            
            best.first = result;
            best.second = _coefs;
            norms = _collapse->function.coefs;
            docopy = false;
          }

          if( docopy )
          {
            _coefs = best.second;
            _collapse->function.coefs = norms;
          }

          // if above were not pre-runs, this is it.
          if( not prerun ) return best.first;

          // Continue runs with a tolerance criterion only.
          types :: t_unsigned itersave = t_Solver :: itermax;
          bool verbsave = t_Solver :: verbose;
          t_Solver :: itermax = 0;
          t_Solver :: verbose = verbose;
          types::t_real result = t_Solver :: operator()( _coefs, _collapse );
          t_Solver :: itermax = itersave;
          t_Solver :: verbose = verbsave;
          return result;
        }

      template< class T_COEFS, class T_TARGETS, class T_COLLAPSE > 
        void bestofcoefs( T_COEFS &_coefs, const T_TARGETS &_targets, T_COLLAPSE &_collapse,
                          types::t_real _howrandom = 5e-1, types::t_unsigned n = 1,
                          bool _verbose = false )
        {
          __DOASSERT( n==0, "0 runs required on input. We are done here.\n" )

          typedef typename std::pair< types::t_real, T_COEFS > t_Best;
          t_Best best;

          bool docopy( false );
          _collapse.update_all( _coefs );
          _collapse.evaluate( _coefs, _targets );
          best.first = _collapse.evaluate( _coefs, _targets );
          best.second = _coefs;
          size_t fieldsize = (size_t) std::ceil( std::log( (types::t_real) n)/std::log(10e0) );
          ++fieldsize;
          if( _verbose ) std::cout << "Try Coefs " << std::setw(fieldsize)
                                   << 1 << ": " << std::setw(15)
                                   << std::fixed
                                   << std::setprecision(3) << best.first << "\n";
          for( types::t_unsigned i(2); i <= n; ++i )
          {
            docopy = true;
            _collapse.create_coefs( _coefs, _howrandom );
            _collapse.update_all( _coefs );
            types::t_real result( _collapse.evaluate( _coefs, _targets ) );

            if( _verbose ) std::cout << "Try Coefs " << std::setw(fieldsize)
                                     << i << ": " << std::setw(15)
                                     << std::fixed
                                     << std::setprecision(3) << result << "\n";

            if( math::geq( result, best.first ) ) continue;
            
            best.first = result;
            best.second = _coefs;
            docopy = false;
          }

          if( docopy ) _coefs = best.second;
        }
  }
}  // namespace LaDa
#endif
