//
//  Version: $Id$
//
#ifndef _SEPARABLES_BESTOF_H_
#define _SEPARABLES_BESTOF_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <utility> 

#include <opt/types.h>
#include <opt/fuzzy.h>
#include <opt/debug.h>


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

    public:
      //! Constructor.
      BestOf ( types::t_unsigned _n = 1 ) : t_Solver(), n(_n), verbose(false) {}
      //! Copy Constructor.
      BestOf ( BestOf &_c ) : t_Solver( _c ), n(_c.n), verbose(_c.verbose) {}
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

        typedef typename std::pair< types::t_real, T_COEFS > t_Best;
        t_Best best;

        bool docopy( false );
        best.first = t_Solver :: operator()( _coefs, _collapse );
        best.second = _coefs;
        if( verbose ) std::cout << "Run 1: " << best.first << "\n";
        for( types::t_unsigned i(1); i <= n; ++i )
        {
          docopy = true;
          _collapse->create_coefs( _coefs );
          types::t_real result( t_Solver ::operator()( _coefs, _collapse ) );

          if( verbose ) std::cout << "Run " << i << ": " << result << "\n";

          if( Fuzzy::geq( result, best.first ) ) continue;
          
          best.first = result;
          best.second = _coefs;
          docopy = false;
        }

        if( docopy ) _coefs = best.second;
        return best.first;
      }
}
#endif
