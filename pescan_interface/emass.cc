//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iterator>
#include <algorithm>
#include <cmath>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <physics/physics.h>
#include <opt/tuple_io.h>

#include "emass.h"

                                                    
namespace LaDa
{
  namespace Pescan
  {
    namespace details
    {
      template< class T_CONTAINER > void double_results( T_CONTAINER& _results );
    }

    template< class T_BANDS >
      void save_bands( const T_BANDS& _bands )
      {
        std::ofstream ofs("_bands");
        boost::archive::text_oarchive oa(ofs);
        oa << _bands;
      }
    template< class T_BANDS >
      void load_bands( T_BANDS& _bands )
      {
        std::ifstream ofs("_bands");
        boost::archive::text_iarchive oa(ofs);
        oa >> _bands;
      }

    void eMass :: operator()
    (
      const Interface& _interface,
      const math::rMatrix3d &_ocell, 
      const Crystal::Structure &_structure,
      const types::t_real &_eref,
      t_Output &_out 
    ) const
    {
      namespace bnu = boost::numeric::ublas;
      __DOASSERT( order >= 2*(npoints-1) + 1,
                  "Interpolation order too large, "
                  "or number of interpolation points to small.\n" )
      __DOASSERT( order < 2,
                  "Interpolation order too small to get second derivative.\n" )
      Interface interface( _interface );

      // setup functional.
      interface.escan.nbstates = nbstates;
      interface.set_method( Interface :: FOLDED_SPECTRUM );
      interface.set_reference( _eref );

      // setup kpoints.
      const bool is_gamma( kpoint.squaredNorm() < 1e-8 );
      const types::t_int nfirst( 1-npoints );
      const types::t_int nlast( is_gamma ? 1: npoints );
      std::vector< math::rVector3d > kpoints;
      for( types::t_int i(nfirst); i < nlast; ++i )
        kpoints.push_back(  kpoint + types::t_real( i ) * direction * stepsize );

      // computes eigenvalues.
      const math::rVector3d ndir( direction.normalized() );
      math::rVector3d korigin( kpoint );
      distort_kpoint( _ocell, _structure.cell, korigin );
      typedef std::pair< types::t_real, std::vector<types::t_real> > t_Eig;
      std::vector<t_Eig> eigs;
      foreach( const math::rVector3d &kp, kpoints )
      {
        compute_
        ( 
          interface,
          distort_kpoint( _ocell, _structure.cell, kp )
        );
        eigs.push_back
        (
          t_Eig
          (
            (interface.escan.kpoint - korigin).dot(ndir),
            interface.eigenvalues
          )
        );
      }
      if( is_gamma ) details::double_results( eigs );
      
      // computes effective masses.
      _out.clear();
      const size_t center_index( eigs.size()>>1 );
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
      for( size_t band(0); band < nbstates; ++band )
      {
        const types::t_real center_eigenvalue( eigs[ center_index ].second[band] );
        // Setups interpolation.
        bnu::matrix<types::t_real> Amatrix( eigs.size(), order + 1);
        bnu::vector<types::t_real> Bvector( eigs.size() );
        bnu::vector<types::t_real> Xvector( order + 1 );
        for( size_t i(0); i < eigs.size(); ++i )
        {
          const types::t_real w( weight_( i, center_index ) );
          Bvector(i) = eigs[i].second[band] * w;
          for( size_t j(0); j < order + 1; ++j )
            Amatrix(i,j) = std::pow( eigs[i].first, int(j) ) * w;
        }

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
        _out.push_back( t_Output::value_type( center_eigenvalue, mass ) );
      } // end of loop over bands.
    } // end of functor. 


    void eMass :: compute_( Interface &_interface, const math::rVector3d &_kpoint ) const
    {
      size_t oldnbstates = _interface.escan.nbstates;
      const bool is_gamma( math::is_zero(_kpoint.squaredNorm()) );
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
    types::t_real eMass :: weight_( size_t _i, size_t _j ) const
    {
      const types::t_int u( std::abs( types::t_int(_i) - types::t_int(_j) ) );
      return u == 0? 1e0: 1e0 / std::pow( types::t_real(u), 4 );
    }

    bool eMass :: load( const TiXmlElement &_node, const std::string &_name )
    {
      const TiXmlElement *parent( opt::find_node( _node, "Functional", "type", _name ) );
      if( not parent ) return false;
      if( parent->Attribute("order") )
        order = boost::lexical_cast< size_t >( parent->Attribute("order") );
      if( parent->Attribute("npoints") )
        npoints = boost::lexical_cast< size_t >( parent->Attribute("npoints") );
      if( parent->Attribute("stepsize") )
        stepsize = boost::lexical_cast< types::t_real >( parent->Attribute("stepsize") );
      if( parent->Attribute("kpoint") )
      {
        boost::tuples::tuple<types::t_real, types::t_real, types::t_real> t;
        opt::tuples::read( parent->Attribute( "kpoint" ), t );
        kpoint[0] = boost::tuples::get<0>( t );
        kpoint[1] = boost::tuples::get<1>( t );
        kpoint[2] = boost::tuples::get<2>( t );
      }
      if( parent->Attribute("direction") )
      {
        boost::tuples::tuple<types::t_real, types::t_real, types::t_real> t;
        opt::tuples::read( parent->Attribute( "direction" ), t );
        direction[0] = boost::tuples::get<0>( t );
        direction[1] = boost::tuples::get<1>( t );
        direction[2] = boost::tuples::get<2>( t );
      }
      if( parent->Attribute("nbstates") )
        nbstates = boost::lexical_cast< size_t >( parent->Attribute("nbstates") );

      cgs.load( *parent );
      return true;
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
