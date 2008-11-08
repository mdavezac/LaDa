//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef _LADADEBUG
# include <boost/bind.hpp>
# include <boost/filesystem/operations.hpp>
#endif
#include <boost/tuple/tuple.hpp>

#include <mpi/mpi_object.h>
#include <opt/initial_path.h>
#include <opt/debug.h>
#include <physics/physics.h>

#include "bandgap.h"
#include "dipole_elements.h"

// declares fortran interface
extern "C"
{
  void FC_FUNC_(iaga_set_cell,IAGA_SET_CELL)
               (const double*, const double[3][3], const int*, const int*, const int* );
  void FC_FUNC_(iaga_dipole_elements, IAGA_DIPOLE_ELEMENTS)
               ( double[], const int*, const int*, const int*,
                 const int*, const char*, const int*, const char*, const int* );
}


namespace LaDa
{
  namespace Pescan
  {
    void dipole_elements( std::vector< Dipole > &_dipoles,
                          const BandGap& _bandgap,
                          const Crystal::Structure &_structure,
                          types::t_real _degeneracy )
    {
      __DODEBUGCODE
      ( 
        namespace bl = boost::lambda; 
        namespace bfs = boost::filesystem; 
      )
      __ASSERT( _bandgap.eigenvalues.size() < 2,
                "Not enough eigenvalues were computed.\n" )
      __ASSERT( _bandgap.bands.gap() < 3e0 * _degeneracy,
                "Band gap too small or degeneracy tolerance too large.\n" )
      __ASSERT
      (
        std::find_if
        ( 
          _bandgap.eigenvalues.begin(), 
          _bandgap.eigenvalues.end(),
          boost::bind( &Fuzzy::gt<types::t_real>, _bandgap.bands.vbm, _1 )
        ) == _bandgap.eigenvalues.end(),
        "VBM and eigenvalues are not coherent.\n"
      )
      __ASSERT
      (
        std::find_if
        ( 
          _bandgap.eigenvalues.begin(), 
          _bandgap.eigenvalues.end(), 
          boost::bind( &Fuzzy::gt<types::t_real>, _bandgap.bands.cbm, _1 )
        ) == _bandgap.eigenvalues.end(),
        "Bandgap and eigenvalues are not coherent.\n"
      )
      
      
      // Finds out the number of conduction and valence bands.
      std::vector<int> bands;
      int ncond(0), nval(0), startval(-1);
      int i = 0;
      foreach( const types::t_real eigenval, _bandgap.eigenvalues )
      {
        if( std::abs( eigenval - _bandgap.bands.vbm ) < _degeneracy )
          bands.push_back( i ), ++nval;
        ++i;
      }
      i = 0;
      foreach( const types::t_real eigenval, _bandgap.eigenvalues )
      {
        if( std::abs( eigenval - _bandgap.bands.cbm ) < _degeneracy )
          bands.push_back( i ), ++ncond;
        ++i;
      }
      
      // Sends cell parameters to fortran
      const int x( boost::tuples::get<0>( _bandgap.genpot.mesh ) );
      const int y( boost::tuples::get<1>( _bandgap.genpot.mesh ) );
      const int z( boost::tuples::get<2>( _bandgap.genpot.mesh ) );
      const types::t_real alat( _structure.scale / Physics::a0("A") );
      FC_FUNC_(iaga_set_cell,IAGA_SET_CELL)
              ( &alat, _structure.cell.x, &x, &y, &z );

      // Call fortran code.
      const size_t dipole_array_size( 2 * 3 * 4 * nval * ncond);
      double dipoles[ dipole_array_size ];
      const boost::filesystem::path holepath
      (
        (
          _bandgap.escan.method == Interface::ALL_ELECTRON ?
            opt::InitialPath::path() / _bandgap.get_dirname():
            opt::InitialPath::path() / _bandgap.get_dirname() / "cbm"
        ) / _bandgap.escan.rspace_wfn
      ); 
      const char* holename = holepath.string().c_str();
      const int holesize( holepath.string().size() );
      const boost::filesystem::path electronpath
      (
        (
          _bandgap.escan.method == Interface::ALL_ELECTRON ?
            opt::InitialPath::path() / _bandgap.get_dirname():
            opt::InitialPath::path() / _bandgap.get_dirname() / "vbm"
        ) / _bandgap.escan.rspace_wfn
      ); 

      __ASSERT( not bfs::exists( holepath ), holepath << " does not exist.\n" )
      __ASSERT( not bfs::exists( electronpath ), electronpath << " does not exist.\n" )
      const char* electronname = electronpath.string().c_str();
      const int electronsize( electronpath.string().size() );
      FC_FUNC_(iaga_dipole_elements, IAGA_DIPOLE_ELEMENTS)
              ( dipoles, &(bands[0]), &(bands[nval]), &nval, &ncond,
                electronname, &electronsize, holename, &holesize );
      std::cout << nval << " " << ncond << " " << bands[0] << " " << bands[nval] << "\n";

      // copies results to output vector.
      // loop positions determined by order of dipoles in fortran code.
      Dipole dipole;
      types::t_real *i_dip( dipoles ); 
      for( size_t l(0); l < nval; ++l )
      {
        dipole.band2band.first = bands[l];
        for( size_t k(0); k < ncond; ++k )
        {
          dipole.band2band.second = bands[nval + k];
          for( size_t sp(0u); sp < 4u; ++sp )
          {
            const types::t_real deltaE(  _bandgap.eigenvalues[ dipole.band2band.second ]
                                        -_bandgap.eigenvalues[ dipole.band2band.first ]  );
            dipole.spin2spin = Dipole::t_SpinTransition( sp );
            for( size_t r(0u); r < 3u; ++r )
            {
              dipole.r[0] = std::complex<types::t_real>( *i_dip, *(i_dip+1) ) * deltaE;
              i_dip += 2;
              dipole.r[1] = std::complex<types::t_real>( *i_dip, *(i_dip+1) ) * deltaE;
              i_dip += 2;
              dipole.r[2] = std::complex<types::t_real>( *i_dip, *(i_dip+1) ) * deltaE;
              i_dip += 2;
            } // loop over cartesian coordinates
            _dipoles.push_back( dipole );
          } // loop over spins
        } // loop over conduction bands
      } // loop over valence bands
    }

    
    types::t_real oscillator_strength( const std::vector<Dipole> &_dipoles )
    {
      types::t_real result;
      foreach( const Dipole &dipole, _dipoles )
        result += std::real(   dipole.r[0] * std::conj(dipole.r[0])
                             + dipole.r[1] * std::conj(dipole.r[1])
                             + dipole.r[2] * std::conj(dipole.r[2]) );
      return result;
    }
    types::t_real oscillator_strength( const BandGap& _bandgap,
                                       const Crystal::Structure &_structure,
                                       types::t_real _degeneracy, bool _print )
    {
      std::vector< Dipole > dipoles;
      dipole_elements( dipoles, _bandgap, _structure, _degeneracy );
      if( _print ) foreach( const Dipole &dipole, dipoles )
                     std::cout << dipole << "\n";
      return oscillator_strength( dipoles );
    }
    std::ostream& operator<<( std::ostream &_stream, const Dipole& _dip )
    {
      const std::string a(    _dip.spin2spin == Dipole::DOWN2DOWN
                           or _dip.spin2spin == Dipole::DOWN2UP ? "-": "+" );
      const std::string b(    _dip.spin2spin == Dipole::DOWN2DOWN
                           or _dip.spin2spin == Dipole::UP2DOWN ? "-": "+" );
      _stream << "<" << a << _dip.band2band.first
              << "| r |" << b << _dip.band2band.second << "> = { "
              << _dip.r[0] << ", " << _dip.r[1] << ", " << _dip.r[2] << " }";
      return _stream;
    }
  }
} // namespace LaDa
