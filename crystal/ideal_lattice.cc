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
      bool operator()( const Eigen::Vector3d& _a, const Eigen::Vector3d &_b )
      {
        const types::t_real a = _a.squaredNorm();
        const types::t_real b = _b.squaredNorm();
        if( std::abs( a - b ) > m ) return a < b;
        if( math::neq( _a[0], _b[0] ) ) return _a[0] < _b[0];
        if( math::neq( _a[1], _b[1] ) ) return _a[1] < _b[1];
        return _a[2] < _b[2];
      }
      types::t_real m;
    };

    void find_first_neighbors( std::vector< Eigen::Vector3d > &_positions,
                               const Eigen::Matrix3d &_cell,
                               const size_t _n )
    {
      __ASSERT( _positions.size() == 0, "Position vector is empty.\n" )
      // first recenters around first atom.
      const Eigen::Matrix3d inv_cell( !_cell );
      const Eigen::Vector3d origin( _positions[0] );
      types::t_real mindist = -1e0;
      foreach( Eigen::Vector3d &pos, _positions )
      {
        pos -= origin;
        pos = inv_cell * pos;
        for( size_t i(0); i < 3; ++i )
          pos[i] = pos[i] - rint(pos[i]);
        pos = _cell * pos;
        const types::t_real d(pos.squaredNorm());
        if( d < types::tolerance ) continue;
        if( mindist > d  or mindist < 0e0 ) mindist = d;
      }
      Eigen::Vector3i range
      (
        std::sqrt( _n * mindist / _cell.col(0).squaredNorm() ),
        std::sqrt( _n * mindist / _cell.col(1).squaredNorm() ),
        std::sqrt( _n * mindist / _cell.col(2).squaredNorm() )
      ); 
      _positions.reserve( _n + (2*range(0)+2) * (2*range(1)+1) * (2*range(2)+1) );
      {
        const std::vector< Eigen::Vector3d > copy( _positions );
        std::vector< Eigen::Vector3d > :: const_iterator i_pos = copy.begin();
        std::vector< Eigen::Vector3d > :: const_iterator i_pos_end = copy.end();
        for(; i_pos != i_pos_end; ++i_pos )
        {
          for( types::t_int i(-range(0) ); i <= range(0); ++i )
            for( types::t_int j(-range(1) ); j <= range(1); ++j )
              for( types::t_int k(-range(2) ); k <= range(2); ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;
                const Eigen::Vector3d d( (*i_pos) + _cell * Eigen::Vector3d( i,j,k )  );
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

    boost::tuples::tuple< Eigen::Matrix3d, Eigen::Vector3d >
      retrieve_deformation( const Structure &_structure, const size_t _nneigs )
      {
        typedef boost::tuples::tuple< Eigen::Matrix3d, Eigen::Vector3d > t_Result;
        __TRYBEGIN
        namespace bt = boost::tuples;
        namespace bnu = boost::numeric::ublas;
        __ASSERT( _structure.lattice == NULL, "No lattice.\n" )
      
        // Computes list of ideal first neighbors.
        std::vector< Eigen::Vector3d > ideals;
        foreach( const Lattice::t_Site &site, _structure.lattice->sites )
          ideals.push_back( site.pos );
        find_first_neighbors( ideals, _structure.lattice->cell, _nneigs );
      
        // Computes list of _structure first neighbors.
        // First looks for barycenter of all atoms.
        // Then picks closest atom of site 0 as the one around which to look
        // for nearest neighbors.
        std::vector< Eigen::Vector3d > non_ideals;
        non_ideals.reserve( _structure.atoms.size() );
        const Eigen::Vector3d origin( _structure.lattice->sites.front().pos );
        Eigen::Vector3d barycenter;
        foreach( const Structure::t_Atom &atom, _structure.atoms )
        {
          non_ideals.push_back( atom.pos );
          barycenter += atom.pos;
        }
        barycenter =  barycenter * 1e0/types::t_real( non_ideals.size() ) ;
        types::t_real mindist(-1);
        types::t_int minindex(-1);
        types::t_int index(0);
        foreach( const Structure::t_Atom &atom, _structure.atoms )
        {
          if( atom.site != 0 ) { ++index; continue; }
          const types::t_real d( (atom.pos - barycenter).squaredNorm() );
          if( mindist > d or minindex < 0 )
          {
            mindist = d;
            minindex = index;
          }
        }
        __DOASSERT( minindex == -1, "Could not idealize structure. Are site indices set?\n" )
        if( minindex != 0 )  std::swap( non_ideals.front(), non_ideals[minindex] );
        std::cout << " ?? " << non_ideals.front() << "\n";
        find_first_neighbors( non_ideals, _structure.cell, _nneigs );
      
        // computes transformation matrix one row at a time. 
        std::cout << "Retrieving ideal lattice.\n";
        Fitting::Cgs cgs;
        cgs.verbose = false;
        cgs.itermax = 100;
        cgs.tolerance = 1e-12;
        t_Result result; 
        for( size_t r(0); r < 3; ++r )
        {
          // Ax = b
          bnu::vector<types::t_real> x(4), b(4);
          bnu::matrix<types::t_real> A(4,4);
          A = bnu::zero_matrix<types::t_real>(4,4);
          // construct vector and matrices.
          for( size_t i(0); i < 4; ++i )
          {
            x(i) = 0e0;
            b(i) = 0e0;
          }
          typedef std::vector< Eigen::Vector3d > :: const_iterator t_cit;
          t_cit i_ideal = ideals.begin();
          t_cit i_ideal_end = ideals.end();
          t_cit i_non_ideal = non_ideals.begin();
          for(; i_ideal != i_ideal_end; ++i_ideal, ++i_non_ideal )
          {
            __ASSERT( i_non_ideal == non_ideals.end(), "Iterator out of range.\n" )
            for( size_t i(0); i < 3; ++i )
            {
              b(i) += (*i_ideal)(r) * (*i_non_ideal)(i);
              for( size_t j(0); j < 3; ++j )
                A(i,j) += (*i_non_ideal)(i) * (*i_non_ideal)(j);
              A(i,3) += (*i_non_ideal)(i);
            }
            for( size_t j(0); j < 3; ++j )
              A(3,j) += (*i_non_ideal)(j);
            A(3,3) += 1e0;
            b(3) += (*i_ideal)(r);
          }
      
          // solves Ax = b
          Fitting::Cgs::t_Return convergence = cgs( A, x, b );
          std::cout << "  _ fit along " << r << ": " << convergence.first << ", "
                    << convergence.second << ", "  << x << "\n"; 
      
          // stores solution.
          for( size_t i(0); i < 3; ++i ) bt::get<0>(result)(r,i) = x(i); 
          bt::get<1>(result)(r) = x(3);
        }
        std::cout << "  deformation matrix:\n" << bt::get<0>( result ) 
                  << "\n  deformation translation: " << bt::get<1>( result ) << "\n\n";
      
        return result;
        __ENDGROUP__
        catch( ... )
        {
          std::cerr << "Error while searching for derformation matrix.\n";
          Eigen::Matrix3d m; m.zero();
          return t_Result( m, Eigen::Vector3d() );
        }
      };

  } // namespace Crystal

} // namespace LaDa
