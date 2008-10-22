//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <mpi/mpi_object.h>
#include <opt/initial_path.h>
#include <opt/debug.h>
#ifdef _LADADEBUG
# include <boost/bind.hpp>
#endif

#include "bandgap.h"
#include "dipole_elements.h"

// declares fortran interface
extern "C"
{
  void FC_FUNC_(iaga_set_cell,IAGA_SET_CELL)
               (const double*, const double*, const int*, const int*, const int* );
  void FC_FUNC_(iaga_dipole_elements, IAGA_DIPOLE_ELEMENTS)
               (double*, int*, int*, int*, const char*, const int* );
}


namespace Pescan
{
  void dipole_elements( std::vector< Dipole > &_dipoles,
                        const BandGap& _bandgap,
                        const Crystal::Structure &_structure,
                        types::t_real _degeneracy )
  {
    __DODEBUGCODE( namespace bl = boost::lambda; )
    int ncond(0), nval(0), startval(-1);
    const size_t dipole_array_size( 2 * 3 * 4 * nval * ncond );
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
    int i = 0;
    foreach( const types::t_real eigenval, _bandgap.eigenvalues )
    {
      if( std::abs( eigenval - _bandgap.bands.vbm ) < _degeneracy )
        bands.push_back( i ), ++nval;
      ++i;
    }
    foreach( const types::t_real eigenval, _bandgap.eigenvalues )
    {
      if( std::abs( eigenval - _bandgap.bands.vbm ) < _degeneracy )
        bands.push_back( i ), ++ncond;
      ++i;
    }
    
    // Sends cell parameters to fortran
    FC_FUNC_(iaga_set_cell,IAGA_SET_CELL)
            ( &_structure.scale, _structure.cell.x,
              &bandgap.genpot.x, &bandgap.genpot.y, &bandgap.genpot.z );
    // Call fortran code.
    double dipoles[ dipole_array_size ];
    const char* filename = bandgap.escan.rspace_wfn.string().c_str();
    const int fsize( bandgap.escan.rspace_wfn.string().size() );
    FC_FUNC_(iaga_dipole_elements, IAGA_DIPOLE_ELEMENTS)
            ( dipoles, &nval, &ncond, filename, &fsize );

    Dipole dipole;
    types::t_real *i_dip( dipoles ); 
    for( size_t sp(0u); sp < 4u; ++sp )
      for( dipole.band2band.first=startval;
           dipole.band2band.first < nval + startval; 
           ++dipole.band2band.first )
        for( dipole.band2band.second = startval + nval;
             dipole.band2band.second < ncond + startval + nval; 
             ++dipole.band2band.second )
          for( size_t r(0u); r < 3u; ++r )
          {
            dipole.spin2spin = Dipole::t_SpinTransition( sp );
            dipole.r[0] = std::complex<types::t_real>( *i_dip, *(i_dip+1) );
            i_dip += 2;
            dipole.r[1] = std::complex<types::t_real>( *i_dip, *(i_dip+1) );
            i_dip += 2;
            dipole.r[2] = std::complex<types::t_real>( *i_dip, *(i_dip+1) );
            _dipoles.push_back( dipole );
          }
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
