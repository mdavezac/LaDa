//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iterator>
#include <algorithm>
#include <cmath>

#include <physics/physics.h>

#include "emass.h"

                                                    
namespace LaDa
{
  namespace Pescan
  {
    namespace details
    {
      template< class T_CONTAINER > void double_results( T_CONTAINER& _results );
    }

    void eMass :: operator()
    (
      const Interface& _interface,
      const atat::rMatrix3d &_ocell, 
      const Crystal::Structure &_structure,
      const atat::rVector3d &_at,
      const atat::rVector3d &_direction,
      const size_t &_nbstates,
      const types::t_real &_eref,
      std::vector< std::pair<types::t_real, types::t_real> > &_out 
    ) const
    {
      namespace bnu = boost::numeric::ublas;
      __DOASSERT( order >= npoints,
                  "Interpolation order too small, "
                  "or number of interpolation points to large.\n" )
      typedef std::vector< std::pair<types::t_real, types::t_real> > t_Results;
      Interface interface( _interface );

      // setup functional.
      interface.escan.nbstates = _nbstates;
      interface.set_method( Interface :: FOLDED_SPECTRUM );
      interface.set_reference( _eref );

      // setup kpoints.
      const bool is_gamma( atat::norm2( _at ) < 1e-8 );
      const types::t_int nfirst( 1-npoints );
      const types::t_int nlast( is_gamma ? 1: npoints );
      std::vector< atat::rVector3d > kpoints;
      for( types::t_int i(nfirst); i < nlast; ++i )
        kpoints.push_back(  _at + types::t_real( i ) * _direction * stepsize );

      // computes eigenvalues.
      const atat::rVector3d ndir( _direction * 1e0 / std::sqrt( atat::norm2( _direction ) ) );
      typedef std::pair< types::t_real, std::vector<types::t_real> > t_Eigs;
      std::vector< t_Eigs > eigs;
      foreach( const atat::rVector3d &kpoint, kpoints )
      {
        t_Eigs eig;
        compute_
        ( 
          interface,
          distort_kpoint( _ocell, _structure.cell, kpoint )
        );
        eigs.push_back
        (
          t_Eigs
          (
            (interface.escan.kpoint - kpoint) * ndir,
            interface.eigenvalues
          )
        );
      }
      if( is_gamma ) details::double_results( eigs );

      
      // computes effective masses.
      _out.clear();
      const size_t center_index( npoints - 1 );
      const types::t_real pi
      (
        3.1415926535897932384626433832795028841971693993751058209749445920e0
      );
      const types::t_real factor
      ( // transforms energy and k vector units to atomic units(Hartree).
        // structure scale is expected in Angstroem.
        // energy from escan should be in eV.
          Physics::Hartree("eV") * 2e0 * pi * pi 
        * Physics::a0("A") * Physics::a0("A") 
        / _structure.scale / _structure.scale 
      );
      for( size_t band(0); band < _nbstates; ++band )
      {
        const types::t_real center_eigenvalue( eigs[ center_index ].second[band] );
        // Setups interpolation.
        bnu::matrix<types::t_real> Amatrix( kpoints.size(), order + 1);
        bnu::vector<types::t_real> Bvector( kpoints.size() );
        bnu::vector<types::t_real> Xvector( order + 1 );
        for( size_t i(0); i < kpoints.size(); ++i )
        {
          Bvector(i) = eigs[i].second[band] * weight_( i, center_index );
          for( size_t j(0); j < order + 1; ++j )
            Amatrix(i,j) = std::pow( eigs[i].first, j ) * weight_(i, center_index);
        }
        for( size_t j(0); j < order + 1; ++j ) Xvector(j) = 0e0;

        // Performs least square fit.
        bnu::matrix<types::t_real> A = bnu::prec_prod( bnu::trans( Amatrix ), Amatrix );
        LaDa::Fitting::Cgs::t_Return
          result = cgs
                   ( 
                     A,
                     Xvector,
                     bnu::prec_prod( bnu::trans( Amatrix ), Bvector )
                   );

        // finally enters result.
        const types::t_real mass( factor / Xvector(2) ); 
        _out.push_back( t_Results::value_type( center_eigenvalue, mass ) );
      } // end of loop over bands.
    } // end of functor. 


    void eMass :: compute_( Interface &_interface, const atat::rVector3d &_kpoint ) const
    {
      size_t oldnbstates = _interface.escan.nbstates;
      const bool is_gamma( atat::norm2( _kpoint ) < 1e-8 );
      if( is_gamma ) _interface.escan.nbstates = oldnbstates >> 1;
      _interface.escan.kpoint = _kpoint;

      _interface();
      if( is_gamma )
      {
        _interface.eigenvalues.reserve( _interface.eigenvalues.size() << 1 );
        std::copy( _interface.eigenvalues.begin(), _interface.eigenvalues.end(),
                   std::back_inserter( _interface.eigenvalues ) );
        _interface.escan.nbstates = oldnbstates;
      }
      std::sort( _interface.eigenvalues.begin(), _interface.eigenvalues.end() );
    }

    //! Adds a weight to interpolation points.
    types::t_real weight_( size_t _i, size_t _j )
    {
      const types::t_int u( std::abs( types::t_int(_i) - types::t_int(_j) ) );
      return 1e0 / types::t_real( std::pow( u, 2 ) );
    }

    namespace details 
    {
      template< class T_CONTAINER > void double_results( T_CONTAINER &_results )
      {
        typedef typename T_CONTAINER :: value_type t_Type;
        typedef typename T_CONTAINER :: const_reverse_iterator t_crit;

        _results.reserve( _results.size() * 2 - 1 );
        t_crit i_first = _results.rbegin();
        t_crit i_end = _results.rend();
        for(++i_first; i_first != i_end; ++i_first )
          _results.push_back( t_Type( -i_first->first, i_first->second ) );
      }
    } // namespace details

  } // namespace Pescan
} // namespace LaDa
