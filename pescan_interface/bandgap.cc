//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <algorithm>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <physics/physics.h>
#include <print/stdout.h>
#include <print/manip.h>
#include <math/fuzzy.h>
#include <opt/debug.h>
#include <opt/initial_path.h>
#include <opt/tinyxml.h>
#include <crystal/lattice.h>

#include "bandgap.h"

namespace LaDa
{
  namespace Pescan
  {

#   ifdef _NOLAUNCH
      //! Fake and fast functional, when compiled with --enable-nolaunch
      void nolaunch_functional( const Crystal::Structure &_str, Bands &bands );
#   endif

    types::t_real BandGap::folded_spectrum(types::t_unsigned _bgstates)
    {
      const t_Path olddirname( dirname );
      create_directory();

      cbm_eigs.clear();
      vbm_eigs.clear();

      // First compute CBM
      computation = CBM;
      escan.Eref = Eref.cbm;
      if ( not do_destroy_dir ) 
      {
        dirname = olddirname / "cbm";
        create_directory();
      }
      create_potential();
      __TRYCODE( launch_pescan();
                 read_and_throw();,
                 "Error while computing CBM.\n") 
      bands.cbm = find_closest_eig( Eref.cbm );
      cbm_eigs = Interface::eigenvalues;

      // Then computes VBM
      computation = VBM;
      escan.Eref = Eref.vbm;
      if ( not do_destroy_dir ) 
      {
        dirname = olddirname / "vbm";
        create_directory();
        create_potential();
      }
      __TRYCODE( launch_pescan();
                 read_and_throw();,
                 "Error while computing VBM.\n") 
      bands.vbm = find_closest_eig( Eref.cbm );
      vbm_eigs = Interface::eigenvalues;

  #ifndef _NOLAUNCH
      if( do_correct and bands.gap() < metallicity ) correct( olddirname );
      if( bands.gap() < metallicity )
      {
        cbm_eigs.clear();
        vbm_eigs.clear();
        dirname = olddirname / "vbm";
        destroy_directory();
        dirname = olddirname / "cbm";
        destroy_directory();
        dirname = olddirname;
        escan.method = ALL_ELECTRON;
        Print :: out << "Failed to find non-metallic band gap.\n"
                     << "Will try an all-electron calculation." << Print::endl;
        all_electron( _bgstates );
      }
  #else
      nolaunch_functional( _str, bands );
  #endif // _NOLAUNCH

      destroy_directory_();
      dirname = olddirname;

      return bands.gap();
    }

    void BandGap :: correct(const t_Path & _dir) 
    {
      Escan saved_state = escan;
      Bands keeprefs = Eref;
      escan.wavefunction_in = escan.wavefunction_out;
      escan.read_in.resize( escan.nbstates );
      std::vector<types::t_unsigned> :: iterator i_r = escan.read_in.begin();
      std::vector<types::t_unsigned> :: iterator i_r_end = escan.read_in.end();
      for(types::t_unsigned u=1; i_r != i_r_end; ++i_r, ++u ) *i_r = u;

      types::t_unsigned i = 0;
      do 
      {
        Print::out __DODEBUGCODE(<< __SPOT_ERROR)
                   << " Found metallic band gap: " << bands.gap() << "\n" 
                   << "    VBM=" << bands.vbm << "  CBM=" << bands.cbm << "\n" 
                   << " Will Try and modify references.\n";

        if( std::abs( bands.vbm - Eref.vbm ) > std::abs( bands.cbm - Eref.cbm ) )
        {
          Eref.vbm -= inc_dec;
          computation = VBM;
          escan.Eref = Eref.vbm;
          Print::out << "Modifying VBM reference: " << Eref.vbm
                     << Print::endl << Print::flush;
          vbm_eigs.clear();
        }
        else
        {
          Eref.cbm += inc_dec;
          computation = CBM;
          escan.Eref = Eref.cbm;
          Print::out << "Modifying CBM reference: " << Eref.cbm
                     << Print::endl << Print::flush;
          cbm_eigs.clear();
        }
        if ( not do_destroy_dir )
          dirname = _dir / ( computation == CBM ?  "/cbm": "/vbm" );
        do_correct = false;
        
        __TRYCODE( launch_pescan();
                   read_and_throw();,
                      "Error while computing"
                   << ( computation == CBM ?  "/cbm": "/vbm" )
                   << "\n" )
        do_correct = true;
        computation == VBM ?  bands.vbm = find_closest_eig( Eref.vbm ):
                              bands.cbm = find_closest_eig( Eref.cbm );
        computation == VBM ?  vbm_eigs = Interface::eigenvalues:
                              cbm_eigs = Interface::eigenvalues;
        Print::out << "After Correction: VBM=" << bands.vbm
                   << " CBM=" << bands.cbm << Print::endl;
        ++i;
        if ( not do_destroy_dir ) dirname = _dir;
      } while ( bands.gap() < metallicity and i < 5 );

      escan = saved_state;
      Eref  = keeprefs;
    }
    
    types::t_real BandGap::all_electron(types::t_unsigned _bgstates)
    {

      types::t_unsigned oldnbstates = escan.nbstates;
      escan.nbstates += _bgstates + 2;
      if (    escan.potential != Escan::SPINORBIT
           or math::is_zero( escan.kpoint.squaredNorm() ) )
       { escan.nbstates >>= 1; _bgstates >>= 1; }
      __TRYCODE( Interface::operator()();,
                 "All-electron calculation failed.\n" )

      escan.nbstates = oldnbstates;

  #ifndef _NOLAUNCH
      foreach( types::t_real ei, eigenvalues )
        std::cout << " " << ei;
      std::cout << " " << _bgstates << "\n";
      bands.cbm = eigenvalues[ _bgstates ];
      bands.vbm = eigenvalues[ _bgstates - 1 ];
      std::cout << "BandGap: " << bands.gap() << " = " 
                << bands.cbm << " - " << bands.vbm << std::endl;
  #else
      LADA_DOASSERT(true, "not implemented")
  #endif // _NOLAUNCH

      destroy_directory_();


      return bands.gap();
    }

    bool BandGap :: Load( const TiXmlElement &_node )
    {
      namespace bfs = boost::filesystem;
      const TiXmlElement *parent
         = opt::find_node( _node, "Functional", "type", "escan" );
      __DOASSERT( not parent, 
                  "Could not find <Functional type=\"escan\"> in input\n"; )

      // Checks for filename attribute.
      bfs::path path;
      TiXmlDocument doc;
      if(  parent->Attribute( "filename" ) )
      {
        path = Print::reformat_home( parent->Attribute( "filename" ) );
        __DOASSERT( not bfs::exists( path ), path.string() + " does not exist.\n" )
        opt::read_xmlfile( path, doc );
        __DOASSERT( not doc.FirstChild( "Job" ),
                    "Root tag <Job> does not exist in " + path.string() + ".\n" )
        parent = opt::find_node( *doc.FirstChildElement( "Job" ),
                                 "Functional", "type", "escan" );
        __DOASSERT( not parent, 
                    "Could not find <Functional type=\"escan\"> in input\n"; )
      }
      
      __TRYCODE( __DOASSERT( not Interface :: Load_( *parent ), "" ),
                 "Encountered error while loading bandgap interface:\n" )
      
      t_method m = escan.method;
      escan.method = ALL_ELECTRON;
      const TiXmlElement *child = parent->FirstChildElement("Reference");
      if( not child ) child = parent->FirstChildElement("References"); 
      if( not child ) return true;
      if(     (not child->Attribute("VBM") )
          and (not child->Attribute("CBM") ) ) return true;
      if(     (not child->Attribute("VBM") )
          or (not child->Attribute("CBM") ) )
      {
        std::cerr __DODEBUGCODE(<< __SPOT_ERROR)
                  << "Found only one VBM/CBM reference when expecting two in input \n"
                  << "Will continue with all electron calculation." << std::endl;
        return true;
      }

      double d;
      child->Attribute( "CBM", &d );
      Eref.cbm = (types::t_real) d;
      child->Attribute( "VBM", &d );
      Eref.vbm = (types::t_real) d;

      if( Eref.gap() < 0.0 )  
      {
        std::cerr __DODEBUGCODE(<< __SPOT_ERROR)
                  << "Found strange VBM/CBM references in input \n"
                  << "VBM = " <<  Eref.vbm << " > CBM = " << Eref.cbm << "\n" 
                  << "Will continue with all electron calculation." << std::endl;
        return true;
      }

      escan.method = m;
      return true;
    }

    void BandGap :: read_and_throw()
    {
      try { read_result(); }
      catch( std::exception &_e )
      { 
        const t_Path filepath( opt::InitialPath::path() / dirname / escan.output.filename() );
        std::ifstream file;
        file.open( filepath.string().c_str(), std::ios_base::in ); 
        
        __DOASSERT( not file.is_open(),
                      "Could not open escan output file "
                   << filepath << "\n" );
        __DOASSERT( file.bad(),
                      "Error while opening escan output file "
                   << filepath << "\n" );
        char cline[256];
        std::string line("");
        std::cout << "\n\n--- ESCAN OUTPUT ---\n";
        while( not file.eof() )
        {
          file.getline( cline, 256 );
          line = cline;
          std::cout << line << "\n";
        }
        std::cout  << "\n\n--- ESCAN OUTPUT ---\n" << std::endl;
        __THROW_ERROR( _e.what() )
      }
    }

#   ifdef _NOLAUNCH
      void nolaunch_functional( const Crystal::Structure &_str, Bands &bands )
      {
        Crystal::Structure::t_kAtoms::const_iterator i_k = _str.k_vecs.begin();
        Crystal::Structure::t_kAtoms::const_iterator i_k_end = _str.k_vecs.end();
        bands.vbm = 0.0; bands.cbm = 0.0;
        bool which = true, sign = true;
        for(; i_k != i_k_end; ++i_k, which = not which )
        {
          types::t_real a0 = std::abs( i_k->type ) * i_k->pos(0);
          types::t_real a1 = std::abs( i_k->type ) * i_k->pos(1);
          types::t_real a2 = std::abs( i_k->type ) * i_k->pos(2);
   
          bands.vbm +=  a0 * a0 / 5.0 - a1 / 15.0 - a2 * a1;
          bands.cbm += 7.0 * cos( a0 * a1 / 7.0 + a2  ); 
        }
   
        bands.cbm += bands.vbm;
   
        if ( math::gt( bands.cbm, bands.vbm ) ) return;
        std::swap( bands.cbm, bands.vbm );
        if ( math::eq(bands.vbm, bands.cbm ) ) bands.cbm += 0.1;
      }
#   endif

  }
} // namespace LaDa
