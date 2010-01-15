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
#include <print/stdout.h>

#include "bandgap.h"
#include "dipole_elements.h"


// declares fortran interface
extern "C"
{
  void FC_FUNC_(momentum, MOMENTUM)
               ( const char*, const int*,
                 const char*, const int*, 
                 const char*, const int*,
                 const int[], const int[],
                 const int*, const int*,
                 double[], const int*, const int* );

}


namespace LaDa
{
  namespace Pescan
  {
    // Computes dipoles, whether all-electron or folded spectrum.
    void compute_dipoles( std::vector< Dipole > &_dipoles,
                          const std::vector<types::t_real> &_vbm_eigs,
                          const std::vector<types::t_real> &_cbm_eigs,
                          const size_t _kramer, const types::t_real _degeneracy,
                          const types::t_real _Evbm, const types::t_real _Ecbm,
                          const std::string& _inputfilename,
                          const boost::filesystem::path &_valence_path,
                          const boost::filesystem::path &_conduction_path );

    void dipole_elements( std::vector< Dipole > &_dipoles,
                          const BandGap& _bandgap,
                          const types::t_real _degeneracy )
    {
      namespace bfs = boost::filesystem; 
      // true if an all electron calculation.
      const std::string filename
      (
           "escan_input."
         + boost::lexical_cast<std::string>(boost::mpi::communicator().rank()) 
      ); 
      const bfs::path path
      (
          opt::InitialPath::path() / _bandgap.get_dirname() 
      ); 
      const size_t kramer( Fuzzy::is_zero( _bandgap.escan.kpoint.squareNorm() ) ? 2: 1 ); 
      if(     _bandgap.escan.method == Interface :: FOLDED_SPECTRUM 
          and bfs::exists( opt::InitialPath::path() /  _bandgap.get_dirname() / "cbm" )  )
      {
        std::cout << (opt::InitialPath::path() /  _bandgap.get_dirname() / "cbm") <<"\n";
        __ASSERT( not bfs::exists( path / "vbm" / filename ), 
                  "Could not find file " + (path / "vbm" / filename).string() )
        __ASSERT( not bfs::exists( path / "cbm" / filename ), 
                  "Could not find file " + (path / "cbm" / filename).string() )
        compute_dipoles
        (
          _dipoles, _bandgap.vbm_eigs, _bandgap.cbm_eigs,
          kramer, _degeneracy, _bandgap.bands.vbm, _bandgap.bands.cbm,
          filename, path / "vbm", path / "cbm"
        );
      }
      else
      {
        __ASSERT( not bfs::exists( path / filename ), 
                  "Could not find file " + (path / filename).string() + ".\n" )
        compute_dipoles
        (
          _dipoles, _bandgap.eigenvalues, _bandgap.eigenvalues,
          kramer, _degeneracy, _bandgap.bands.vbm, _bandgap.bands.cbm,
          filename, path, path
        );
      }
    }

    void compute_dipoles( std::vector< Dipole > &_dipoles,
                          const std::vector<types::t_real> &_vbm_eigs,
                          const std::vector<types::t_real> &_cbm_eigs,
                          const size_t _kramer,
                          const types::t_real _degeneracy,
                          const types::t_real _Evbm, const types::t_real _Ecbm,
                          const  std::string &_inputfilename,
                          const  boost::filesystem::path &_valence_path,
                          const  boost::filesystem::path &_conduction_path )
    {
      namespace bfs = boost::filesystem; 
      __ASSERT( not bfs::exists( _valence_path ), _valence_path << " does not exist.\n" )
      __ASSERT( not bfs::exists( _conduction_path ), _conduction_path << " does not exist.\n" )
      
      // Finds out the number of conduction and valence bands.
      std::vector<int> valence_indices;
      std::vector<types::t_real> :: const_iterator i_band = _vbm_eigs.begin();
      std::vector<types::t_real> :: const_iterator i_band_end = _vbm_eigs.end();
      for( int index(1) ; i_band != i_band_end; ++i_band, ++index )
        if( std::abs( *i_band - _Evbm ) < _degeneracy  )
          valence_indices.push_back( index );
      __ASSERT( valence_indices.size() == 0,  "VBM and eigenvalues are not coherent.\n" )
      std::vector<int> conduction_indices;
      i_band = _cbm_eigs.begin();
      i_band_end = _cbm_eigs.end();
      for( int index(1) ; i_band != i_band_end; ++i_band, ++index )
        if( std::abs( *i_band - _Ecbm ) < _degeneracy  )
          conduction_indices.push_back( index );
      __ASSERT( conduction_indices.size() == 0,  "CBM and eigenvalues are not coherent.\n" )
      
      // Call fortran code.
      const int ninputpath( _inputfilename.size() );
      const char* char_inputpath( _inputfilename.c_str() );
      const int nwfn_vbm( _valence_path.string().size() );
      const char* char_wfnvbm( _valence_path.string().c_str() );
      const int nwfn_cbm( _conduction_path.string().size() );
      const char* char_wfncbm( _conduction_path.string().c_str() );
      const int nval( valence_indices.size() );
      const int ncond( conduction_indices.size() );
      const int ndip2( nval * _kramer );
      const int ndip3( ncond * _kramer );
      double dipoles[ size_t(ndip3) ][ size_t(ndip2) ][3][2];
      
      FC_FUNC_(momentum, MOMENTUM)
              ( char_inputpath, &ninputpath,
                char_wfnvbm, &nwfn_vbm, 
                char_wfncbm, &nwfn_cbm,  
                &valence_indices[0], &conduction_indices[0],
                &nval, &ncond, &(dipoles[0][0][0][0]), &ndip2, &ndip3 );

      // copies results to output vector.
      // loop positions determined by order of dipoles in fortran code.
      Dipole dipole;
      for( int k(0); k < ndip2; ++k )
      {
        dipole.band2band.first = valence_indices[k / _kramer ] - 1;
        for( int l(0); l < ndip3; ++l )
        {
          dipole.band2band.second = conduction_indices[l / _kramer ] - 1;

          for( size_t r(0u); r < 3u; ++r )
          {
            dipole.p[r] = std::complex<types::t_real>
                          ( 
                            dipoles[l][k][r][0],
                            dipoles[l][k][r][1]
                          );
          } // loop over cartesian coordinates
          _dipoles.push_back( dipole );
        } // loop over conduction bands
      } // loop over valence bands
    }
    
    types::t_real oscillator_strength( const std::vector<Dipole> &_dipoles )
    {
      types::t_real result(0);
      foreach( const Dipole &dipole, _dipoles )
        result += std::real(   dipole.p[0] * std::conj(dipole.p[0])
                             + dipole.p[1] * std::conj(dipole.p[1])
                             + dipole.p[2] * std::conj(dipole.p[2]) );
      return result;
    }
    types::t_real oscillator_strength( const BandGap& _bandgap,
                                       types::t_real _degeneracy, bool _print )
    {
      std::vector< Dipole > dipoles;
      dipole_elements( dipoles, _bandgap, _degeneracy );
      _print = true;
      const types::t_real result( oscillator_strength( dipoles ) );
      if( _print )
      {
        foreach( const Dipole &dipole, dipoles )
          std::cout << dipole << "\n";
        std :: cout << "oscillator strength: " << result << Print::endl;
      }
      return result;
    }
    std::ostream& operator<<( std::ostream &_stream, const Dipole& _dip )
    {
      _stream << "<" << _dip.band2band.first
              << "| p |" << _dip.band2band.second << "> = { "
              << _dip.p[0] << ", " << _dip.p[1] << ", " << _dip.p[2] << " }";
      return _stream;
    }


  } // namespace Pescan
} // namespace LaDa
