//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <minimizer/cgs.h>
#include <opt/debug.h>
#include <opt/fuzzy.h>

#include "ideal_lattice.h"

namespace LaDa
{

  namespace Crystal 
  {

    //! \brief Strict weak ordering for distance from origin.
    //! \details Comparison with origin is always false.
    //!          This way, ends up a the end when sorting.
    struct Order
    {
      Order( types::t_real r ) : m(r) {}
      Order( const Order &_c ) : m(_c.m) {}
      bool operator()( const atat::rVector3d& _a, const atat::rVector3d &_b )
      {
        const types::t_real a = atat::norm2( _a );
        const types::t_real b = atat::norm2( _b );
        if( std::abs( a - b ) > m ) return a < b;
        if( Fuzzy::neq( _a[0], _b[0] ) ) return _a[0] < _b[0];
        if( Fuzzy::neq( _a[1], _b[1] ) ) return _a[1] < _b[1];
        return _a[2] < _b[2];
      }
      types::t_real m;
    };

    void find_first_neighbors( std::vector< atat::rVector3d > &_positions,
                               const atat::rMatrix3d &_cell,
                               const size_t _n )
    {
      __ASSERT( _positions.size() == 0, "Position vector is empty.\n" )
      // first recenters around first atom.
      const atat::rMatrix3d inv_cell( !_cell );
      const atat::rVector3d origin( _positions[0] );
      types::t_real mindist = -1e0;
      foreach( atat::rVector3d &pos, _positions )
      {
        pos -= origin;
        pos = inv_cell * pos;
        for( size_t i(0); i < 3; ++i )
          pos[i] = pos[i] - rint(pos[i]);
        pos = _cell * pos;
        const types::t_real d( atat::norm2( pos ) );
        if( d < types::tolerance ) continue;
        if( mindist > d  or mindist < 0e0 ) mindist = d;
      }
      atat::iVector3d range
      (
        std::sqrt( _n * mindist / atat::norm2( _cell.get_column(0) ) ),
        std::sqrt( _n * mindist / atat::norm2( _cell.get_column(1) ) ),
        std::sqrt( _n * mindist / atat::norm2( _cell.get_column(2) ) )
      ); 
      _positions.reserve( _n + (2*range(0)+2) * (2*range(1)+1) * (2*range(2)+1) );
      {
        const std::vector< atat::rVector3d > copy( _positions );
        std::vector< atat::rVector3d > :: const_iterator i_pos = copy.begin();
        std::vector< atat::rVector3d > :: const_iterator i_pos_end = copy.end();
        for(; i_pos != i_pos_end; ++i_pos )
        {
          for( types::t_int i(-range(0) ); i <= range(0); ++i )
            for( types::t_int j(-range(1) ); j <= range(1); ++j )
              for( types::t_int k(-range(2) ); k <= range(2); ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;
                const atat::rVector3d d( (*i_pos) + _cell * atat::rVector3d( i,j,k )  );
                _positions.push_back( d );
              }
        }
      }

      // then resizes list of atoms to _n 
      std::partial_sort( _positions.begin(), _positions.begin() + _n + 1, 
                         _positions.end(), Order( mindist * 0.25 ) );
      std::swap( _positions.front(), _positions[_n] );
      _positions.resize( _n );
    }

    atat::rMatrix3d retrieve_deformation( const Structure &_structure,
                                          const size_t _nneigs )
    {
      __TRYBEGIN
      namespace bnu = boost::numeric::ublas;
      __ASSERT( _structure.lattice == NULL, "No lattice.\n" )

      // Computes list of ideal first neighbors.
      std::vector< atat::rVector3d > ideals;
      foreach( const Lattice::t_Site &site, _structure.lattice->sites )
        ideals.push_back( site.pos );
      find_first_neighbors( ideals, _structure.lattice->cell, _nneigs );
      std::cout << "ideal: \n";
      foreach( const atat::rVector3d &a, ideals )
        std::cout << " " << a << "\n";
      std::cout << "\n";

      // Computes list of _structure first neighbors.
      std::vector< atat::rVector3d > non_ideals;
      non_ideals.reserve( _structure.atoms.size() );
      const atat::rVector3d origin( _structure.lattice->sites.front().pos );
      types::t_real mindist(-1);
      types::t_int minindex(-1);
      types::t_int index(0);
      foreach( const Structure::t_Atom &atom, _structure.atoms )
      {
        non_ideals.push_back( atom.pos );
        const types::t_real d( atat::norm2( atom.pos - origin ) );
        if( mindist > d or minindex < 0 )
        {
          mindist = d;
          minindex = index;
        }
        ++index;
      }
      __DOASSERT
      (
        mindist > types::tolerance, 
        "Cannot work with translated structure.\n" 
      )
      if( minindex != 0 )  std::swap( non_ideals.front(), non_ideals[minindex] );
      find_first_neighbors( non_ideals, _structure.cell, _nneigs );
      std::cout << "non ideal: \n";
      foreach( const atat::rVector3d &a, non_ideals )
        std::cout << " " << a << "\n";
      std::cout << "\n";

      // computes transformation matrix one row at a time. 
      Fitting::Cgs cgs;
      cgs.verbose = false;
      cgs.itermax = 100;
      cgs.tolerance = 1e-12;
      atat::rMatrix3d result; 
      for( size_t r(0); r < 3; ++r )
      {
        // Ax = b
        bnu::vector<types::t_real> x(3), b(3);
        bnu::matrix<types::t_real> A(3,3);
        A = bnu::zero_matrix<types::t_real>(3,3);
        // construct vector and matrices.
        for( size_t i(0); i < 3; ++i )
        {
          x(i) = 0e0;
          b(i) = 0e0;
        }
        std::vector< atat::rVector3d > :: const_iterator i_ideal = ideals.begin();
        std::vector< atat::rVector3d > :: const_iterator i_ideal_end = ideals.end();
        std::vector< atat::rVector3d > :: const_iterator i_non_ideal = non_ideals.begin();
        for(; i_ideal != i_ideal_end; ++i_ideal, ++i_non_ideal )
        {
          __ASSERT( i_non_ideal == non_ideals.end(), "Iterator out of range.\n" )
          for( size_t i(0); i < 3; ++i )
          {
            b(i) += (*i_ideal)(r) * (*i_non_ideal)(i);
            for( size_t j(0); j < 3; ++j ) A(i,j) += (*i_non_ideal)(i) * (*i_non_ideal)(j);
          }
        }

        // solves Ax = b
        Fitting::Cgs::t_Return convergence = cgs( A, x, b );
        std::cout << r << ": " << convergence.first << ", "
                  << convergence.second << ", "  << x << "\n"; 

        // stores solution.
        for( size_t i(0); i < 3; ++i ) result(r,i) = x(i); 
      }

      return result;
      __ENDGROUP__
      catch( ... )
      {
        std::cerr << "Error while searching for derformation matrix.\n";
        atat::rMatrix3d result; 
        result.zero();
        return result;
      }
    };

  } // namespace Crystal

} // namespace LaDa
