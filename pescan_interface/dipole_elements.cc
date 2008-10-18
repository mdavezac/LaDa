//
//  Version: $Id: interface.h 816 2008-10-17 01:29:20Z davezac $
//

# ifdef HAVE_CONFIG_H
#   include <config.h>
# endif

#include<opt/physics.h>

#include"dipole_element.h"

namespace Pescan
{
  types::t_real dipole_elements( const BandGap& _bandgap,
                                 const Crystal::Structure &_str,
                                 types::t_real _degeneracy )
  {
    if( _bandgap.comm().rank() ) return 0;

    const boost::path filepath(   opt::InitialPath::path()
                                / _bandgap.get_dirname() / "mxmat.d" );
    __ASSERT( _bandgap.eigenvalues.size() < 2,
              "Not enough eigenvalues were computed.\n" )
    __ASSER( bands.gap() > 3e0 * _degeneracy,
             "Band gap too small or degeneracy tolerance too large.\n" )
    __ASSERT
    (
      std::find
      ( 
        _bandgap.eigenvalues.begin(), 
        _bandgap.eigenvalues.end(),
        bl::bind
        ( 
          &Fuzzy::gt<types::t_real>,
          bl::constant( _bandgap.bands.bbm ),
          bl::_1
        ) 
      ) == _bandgap.eigenvalues.end(),
      "VBM and eigenvalues are not coherent.\n"
    )
    __ASSERT
    (
      std::find
      ( 
        _bandgap.eigenvalues.begin(), 
        _bandgap.eigenvalues.end(),
        bl::bind
        ( 
          &Fuzzy::gt<types::t_real>,
          bl::constant( _bandgap.bands.cbm ),
          bl::_1
        ) 
      ) == _bandgap.eigenvalues.end(),
      "Bandgap and eigenvalues are not coherent.\n"
    )


    // Finds out the number of conduction and valence bands.
    size_t ncond(0), nval(0);
    foreach( const types::t_real eigenval, _bandgap.eigenvalues )
      if( std::abs( eigenval - bands.vbm ) < _degeneracy ) ++nval;
      else if( std::abs( eigenval - bands.cbm ) < _degeneracy ) ++ncond;

    // Now writes the dipole moment file,
    std::ofstream file( filepath.string().c_str(), std::ios_base::out|std::ios_base::trunc );

    file << structure.scale / Physics::a0("A") << "\n";
    // prints cell vectors in units of a0 and other
    // whatever nanopes other may be
    for( types::t_unsigned i = 0; i < 3; ++i )
      file << std::fixed << std::setprecision(7) 
           << std::setw(12) << std::setprecision(8)
           << structure.cell(0,i) * structure0.scale / Physics::a0("A")
           << std::setw(12) << std::setprecision(8)
           << structure.cell(1,i) * structure0.scale / Physics::a0("A")
           << std::setw(12) << std::setprecision(8)
           << structure.cell(2,i) * structure0.scale / Physics::a0("A")
           << std::setw(18) << std::setprecision(8) << structure.cell(0,i) 
           << std::setw(12) << std::setprecision(8) << structure.cell(1,i) 
           << std::setw(12) << std::setprecision(8) << structure.cell(2,i) << "\n";
    // prints number of valence bands.
    file << nval << "\n";


    file.flush();
    file.close();
  }
}
