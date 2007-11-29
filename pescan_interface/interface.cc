//
//  Version: $Id$
//
#include "interface.h"
#include <mpi/mpi_object.h>
#include <sstream>

#include <physics/physics.h>
#include <print/manip.h>
#include <print/stdout.h>

#ifdef _NOLAUNCH
#include<eo/utils/eoRNG.h>
#endif

namespace Pescan
{
  types::t_real Interface :: operator()( Ising_CE::Structure &_str)
  {
    if ( escan.method == Escan::FOLDED_SPECTRUM )
      computation = CBM;
      
    std::string olddirname = dirname;
    escan.scale = _str.scale;
    if ( escan.method == Escan::FOLDED_SPECTRUM and (not do_destroy_dir) )
    {
      dirname = olddirname;
      create_directory();
      dirname = olddirname + ( computation == CBM ?  "/cbm": "/vbm" );
    }
    create_directory();
    create_potential();
    launch_pescan( _str );
    types::t_real result = read_result( _str );

    if ( escan.method == Escan::ALL_ELECTRON )
    {
      destroy_directory_();
      escan.Eref = bands;
      return result;
    }
    bands.cbm = result;

    computation = VBM;
    destroy_directory_();
    if ( escan.method == Escan::FOLDED_SPECTRUM and (not do_destroy_dir) )
      dirname = olddirname + ( computation == CBM ?  "/cbm": "/vbm" );
    create_directory();
    create_potential();
    launch_pescan( _str );
    bands.vbm = read_result( _str );
    destroy_directory_();


    if ( bands.gap() > 0.001 )  return bands.gap();

    Bands keeprefs = escan.Eref;
    do
    {
      Print::out << " Found metallic band gap!! " << result << "\n" 
                 << " Will Try and modify references. \n";

      if( std::abs( bands.vbm - escan.Eref.vbm ) > std::abs( bands.cbm - escan.Eref.cbm ) )
      {
        escan.Eref.vbm -= 0.05;
        computation = VBM;
      }
      else
      {
        escan.Eref.cbm += 0.05;
        computation = CBM;
      }
      if ( escan.method == Escan::FOLDED_SPECTRUM and (not do_destroy_dir) )
        dirname = olddirname + ( computation == CBM ?  "/cbm": "/vbm" );
      create_directory();
      create_potential();
      launch_pescan( _str );
      computation == VBM ?  bands.vbm = read_result( _str ):
                            bands.cbm = read_result( _str );
    }
    while ( bands.gap() < 0.001 );

    return bands.gap();
  }
  void Interface :: create_directory()
  {
    std::ostringstream sstr;
    sstr << "mkdir " << dirname;
    system( sstr.str().c_str() );
  }
  void Interface :: destroy_directory_()
  {
    std::ostringstream sstr;
    sstr << "rm -rf " << dirname;
    system( sstr.str().c_str() );
  }
  void Interface :: create_potential()
  {
    write_genpot_input();

    std::ostringstream sstr;
    sstr << "cp " << atom_input << " ";
#ifndef _NOLAUNCH
    sstr <<  genpot.launch << " ";
#endif
    std::vector<std::string> :: const_iterator i_str = genpot.pseudos.begin();
    std::vector<std::string> :: const_iterator i_str_end = genpot.pseudos.end();
    for(; i_str != i_str_end; ++i_str)
      sstr << *i_str << " ";
    sstr << " " << dirname;
    system( sstr.str().c_str() );

    
    sstr.str("");
    sstr << "cd " << dirname << "; ./" << Print::StripDir(genpot.launch);
#ifndef _NOLAUNCH
    system(sstr.str().c_str());
#endif
  }
  types::t_real Interface :: launch_pescan( Ising_CE::Structure &_str )
  {
    std::ostringstream sstr;
    write_escan_input( _str );

    sstr << "cp maskr "; 
#ifndef _NOLAUNCH
    sstr << escan.launch << " ";
#endif
    std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
    std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
    std::vector<std::string> alreadythere;
    for( ; i_so != i_so_end; ++i_so)
    {
      if(    std::find( alreadythere.begin(), alreadythere.end(), i_so->filename ) 
          == alreadythere.end() )
      {
        sstr << i_so->filename << " ";
        alreadythere.push_back(i_so->filename);
      }
    }
    sstr << dirname;
    system( sstr.str().c_str() );

    std::string output = Print::StripEdges(escan.output);
    sstr.str("");
    sstr << "cd " << dirname << "; ./" << Print::StripDir(escan.launch) << " > " << output;
#ifndef _NOLAUNCH
    system(sstr.str().c_str());
#endif
    return 0.0;
  }


  bool Interface :: Load (const TiXmlElement &_node )
  {
    const TiXmlElement *child, *parent;
    std::string str;

    // This whole section tries to find a <Functional type="escan"> tag
    // in _element or its child
    str = _node.Value();
    if ( str.compare("Functional" ) != 0 )
      parent = _node.FirstChildElement("Functional");
    else
      parent = &_node;

    
    while (parent)
    {
      str = "";
      if ( parent->Attribute( "type" )  )
        str = parent->Attribute("type");
      if ( str.compare("escan" ) == 0 )
        break;
      parent = parent->NextSiblingElement("Functional");
    }
    if ( not parent )
    {
      std::cerr << "Could not find an <Functional type=\"escan\"> tag in input file" 
                << std::endl;
      return false;
    } 

    // lester needs to load "module fftw" when in interactive multi-proc mode
    if ( parent->Attribute("system") )
    {
      std::string str = parent->Attribute("system");
      system( str.c_str() );
    }

    do_destroy_dir = true;
    if( parent->Attribute("keepdirs") ) do_destroy_dir = false;

    child = parent->FirstChildElement("GenPot");
    if (    (not child)
         or (not child->Attribute("x")) 
         or (not child->Attribute("y")) 
         or (not child->Attribute("z")) 
         or (not child->Attribute("cutoff")) 
         or (not child->FirstChildElement("Pseudo")) )
    {
      std::cerr << "Need genpot input !!" << std::endl;
      return false;
    }
    child->Attribute("x", &genpot.x);
    child->Attribute("y", &genpot.y);
    child->Attribute("z", &genpot.z);
    child->Attribute("cutoff", &genpot.cutoff);
    if ( child->Attribute("launch") )
      genpot.launch = child->Attribute("launch");
    child = child->FirstChildElement("Pseudo");
    for(; child; child = child->NextSiblingElement() )
      if ( child->Attribute("filename") )
        genpot.pseudos.push_back( child->Attribute("filename") );
    if( not genpot.check() )
    {
      std::cerr << "Insufficient or Incorrect Genpot input!! " << std::endl;
      return false;
    }


    child = parent->FirstChildElement("Wavefunctions");
    if( child and child->Attribute("in") )
      escan.wavefunction_in = child->Attribute("in");
    if( child and child->Attribute("out") )
      escan.wavefunction_out = child->Attribute("out");
    types::t_int j;
    if( parent->Attribute("method", &j) )
      escan.method = ( j == 1 ) ? Escan::FOLDED_SPECTRUM: Escan::ALL_ELECTRON;
    child = parent->FirstChildElement("References");
    if( child )
    {
      if( child->Attribute("VBM") )
        child->Attribute("VBM", &escan.Eref.vbm);
      if( child->Attribute("CBM") )
        child->Attribute("CBM", &escan.Eref.cbm);
    }

    child = parent->FirstChildElement("Hamiltonian");
    if (    ( not child )
         or (not child->Attribute("kinscal"))
         or (not child->Attribute("realcutoff")) 
         or (not child->Attribute("potential")) )
    {
      std::cerr << "Please Specify hamiltonian on input" << std::endl;
      exit(0);
    }
    if ( child->Attribute("launch") )
      escan.launch = child->Attribute("launch");
    child->Attribute("kinscal", &escan.kinscal);
    if( child->Attribute("smooth") )
      child->Attribute("smooth", &escan.smooth);
    if( child->Attribute("potential") )
    {
      child->Attribute("potential", &j);
      switch( j )
      {
        case Escan::NOPOT: std::cerr << "Error, incorrect escan potential " << std::cerr; return false;
        case Escan::LOCAL: escan.potential = Escan::LOCAL;
        case Escan::NONLOCAL: escan.potential = Escan::NONLOCAL;
        case Escan::SPINORBIT: escan.potential = Escan::SPINORBIT;
      }
    }
    if( child->Attribute("realcutoff") )
      child->Attribute("realcutoff", &escan.rcut);
    if( child->Attribute("nbstates") )
      child->Attribute("nbstates", &escan.nbstates);
    child = child->FirstChildElement("SpinOrbit");
    for(; child; child = child->NextSiblingElement("SpinOrbit") )
      if ( child->Attribute("filename") and child->Attribute("izz") )
      {
        SpinOrbit so;
        so.filename = child->Attribute("filename");
        so.izz = child->Attribute("izz");
        if( child->Attribute("s") )
          child->Attribute("s", &so.s);
        if( child->Attribute("p") )
          child->Attribute("p", &so.p);
        if( child->Attribute("d") )
          child->Attribute("d", &so.d);
        if( child->Attribute("pnl") )
          child->Attribute("pnl", &so.pnl);
        if( child->Attribute("dnl") )
          child->Attribute("pnl", &so.dnl);
        escan.spinorbit.push_back(so);
      }
    child = parent->FirstChildElement("Minimizer");
    if( not child )
     return true;

    if( child->Attribute("niter") )
      child->Attribute("niter", &escan.niter);
    if( child->Attribute("nlines") )
      child->Attribute("nlines", &escan.nlines);
    if( child->Attribute("tolerance") )
      child->Attribute("tolerance", &escan.tolerance);

    return true;
  }

  void Interface::write_genpot_input()
  {
    std::ofstream file;
    std::string name = Print::StripEdges(dirname) + "/" + Print::StripEdges(genpot.filename);
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 

    file << Print::StripDir(atom_input) << std::endl 
         << genpot.x << " " 
         << genpot.y << " " 
         << genpot.z << std::endl 
         << genpot.x << " " 
         << genpot.y << " " 
         << genpot.z << std::endl << " 0 0 0 " << std::endl
         << genpot.cutoff << std::endl
         << genpot.pseudos.size() << std::endl;
    std::vector< std::string > :: const_iterator i_str = genpot.pseudos.begin();
    std::vector< std::string > :: const_iterator i_str_end = genpot.pseudos.end();
    for(; i_str != i_str_end; ++i_str )
     file << Print::StripDir(*i_str) << std::endl;
    file.flush();
    file.close();
  }

  void Interface::write_escan_input( Ising_CE::Structure &_str ) 
  {
    types::t_unsigned nbwfn = escan.nbstates; 
    std::ofstream file;
    std::string name = Print::StripEdges(dirname) + "/" + Print::StripEdges(escan.filename);
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    file << "1 " << Print::StripDir(dirname, genpot.output) << "\n"
         << "2 " << wfn_name( escan.wavefunction_out ) << "\n"
         << "3 " << escan.method << "\n";
    switch (computation)
    { 
      case VBM: file << "4 " << escan.Eref.vbm << " "; break;
      case CBM: file << "4 " << escan.Eref.cbm << " "; break;
    }
    file << genpot.cutoff << " " 
         << escan.smooth << " "
         << escan.kinscal << "\n";

    if( escan.method == Escan::ALL_ELECTRON )
    {
      Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
      Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
      for(nbwfn=0; i_atom != i_atom_end; ++i_atom)
      {
        Ising_CE::StrAtom atom; 
        _str.lattice->convert_Atom_to_StrAtom( *i_atom, atom );
        nbwfn += Physics::Atomic::Charge( atom.type );
      }
      if (    escan.potential != Escan::SPINORBIT
           or atat::norm2(escan.kpoint) < types::tolerance ) nbwfn /= 2;
      ++nbwfn; // Adds 1 wavefunction to include CBM
    }
    file << "5 " << nbwfn << "\n"
         << "6 " << escan.niter << " " << escan.nlines << " " << escan.tolerance << "\n";
    
    // Input wavefunctions
    if( do_input_wavefunctions )
    {
      if( escan.method == FOLDED_SPECTRUM )
      {
        file << "7 " << nbwfn << "\n8 ";
        for( types::t_unsigned u=1; u < nbwfn; ++u )
          file << u << ", ";
        file << nbwfn << "\n9 " << wfn_name( escan.wavefunction_in ) << "\n";
      }
      else file << "7 2\n8 " << nbwfn-1 << ", "<< nbwfn
                << "\n9 " << wfn_name( escan.wavefunction_in ) << "\n";
    }
    else  file << "7 0\n8 0\n9 dummy\n10 0 4 4 4 graph.R\n";

    if ( atat::norm2( escan.kpoint ) < types::tolerance )
     file << "11 0 0 0 0 ";
    else 
     file << "11 1 " << escan.kpoint[0] << " " << escan.kpoint[1] << " " << escan.kpoint[2] << " ";
    file << escan.scale << "\n"
         << "12 " << escan.potential << "\n"
         << "13 " << atom_input << "\n"
         << "14 " << escan.rcut << "\n"
         << "15 " << escan.spinorbit.size() << "\n";
    std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
    std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
    for(types::t_unsigned i=16; i_so != i_so_end; ++i_so, ++i )
    {
      file << i << " " << i_so->filename << " " 
           << i_so->izz << " " 
           << i_so->s << " " << i_so->p << " " << i_so->d << " " 
           << i_so->pnl << " " << i_so->dnl << "\n";
    }
    file.flush();
    file.close();
  }

  types::t_real Interface :: read_result( Ising_CE::Structure &_str )
  {
#ifdef _NOLAUNCH
    nolaunch_functional( _str, bands );
    if ( escan.method == Escan :: ALL_ELECTRON ) 
      return bands.gap();
    if ( computation == CBM )
      return bands.cbm;
    return bands.vbm;
#endif
    std::ifstream file;
    std::ostringstream sstr;
    sstr << dirname; 
    sstr << "/" << Print::StripEdges(escan.output);
    file.open( sstr.str().c_str(), std::ios_base::in ); 
    if( file.fail() )
    {
      std::string filename = sstr.str();
      sstr.str("");
      sstr << "Could not open file "
           << filename
           << " in Pescan::Interface::read_result " << std::endl;
      throw std::runtime_error( sstr.str() );
    }

    char cline[256];
    std::string line("");
    while (     ( not file.eof() )
            and line.find("FINAL eigen") == std::string::npos ) 
    {
      file.getline( cline, 256 );
      line = cline;
    }

    types :: t_real eigenvalue;
    if (       escan.nbstates == 1  
          and  escan.method == Escan::FOLDED_SPECTRUM )
    {
      file >> eigenvalue;
      return eigenvalue;
    }

    std::vector<types::t_real> eigenvalues;
    while(     ( not file.eof() )
           and (file >> eigenvalue).good() )
      eigenvalues.push_back( eigenvalue );

    if( escan.method == Escan::ALL_ELECTRON )
    {
      types::t_unsigned n = 0;
      Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
      Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom)
      {
        Ising_CE::StrAtom atom; 
        _str.lattice->convert_Atom_to_StrAtom( *i_atom, atom );
        n += Physics::Atomic::Charge( atom.type );
      }
      if (    escan.potential != Escan::SPINORBIT
           or atat::norm2(escan.kpoint) < types::tolerance )
        n /= 2;
      // watch out! C arrays start at index = 0
      if ( eigenvalues.size() < n or n < 2 )
      {
        std::cerr << "Error, not enough states were computed to determine band gap "
                  << eigenvalues.size() << " vs " << n << " required " << std::endl;
        return -1.0;
      }
      bands.vbm = eigenvalues[n-2];
      bands.cbm = eigenvalues[n-1];
      return bands.gap();
    }

    types :: t_real reference = 0;
    switch (computation)
    { 
      case VBM: reference = escan.Eref.vbm; break;
      case CBM: reference = escan.Eref.cbm; break;
    };
    std::vector<types::t_real> :: const_iterator i_eig = eigenvalues.begin();
    std::vector<types::t_real> :: const_iterator i_eig_end = eigenvalues.end();
    std::vector<types::t_real> :: const_iterator i_eig_result = i_eig;
    types::t_real mini = std::abs(*i_eig-reference); 
    for(++i_eig; i_eig != i_eig_end; ++i_eig)
      if( std::abs( *i_eig - reference ) < mini )
      {
        mini = std::abs( *i_eig - reference );
        i_eig_result = i_eig;
      }
    return *i_eig_result;
  }

#ifdef _NOLAUNCH
  void nolaunch_functional( const Ising_CE::Structure &_str, Interface::Bands &bands )
  {
    Ising_CE::Structure::t_kAtoms::const_iterator i_k = _str.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms::const_iterator i_k_end = _str.k_vecs.end();
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

    if ( bands.vbm > bands.cbm ) bands.swap();
    if ( std::abs( bands.vbm - bands.cbm ) < types::tolerance )
      bands.cbm += 0.1;

  }
#endif
 
}

#ifdef _MPI

namespace mpi
{
  template<>
  bool BroadCast :: serialize<Pescan::Interface::SpinOrbit>( Pescan::Interface::SpinOrbit &_sp )
  {
    if( not serialize( _sp.filename ) ) return false;
    if( not serialize( _sp.izz ) ) return false;
    if( not serialize( _sp.s ) ) return false;
    if( not serialize( _sp.p ) ) return false;
    if( not serialize( _sp.d ) ) return false;
    if( not serialize( _sp.pnl ) ) return false;
    return serialize( _sp.dnl );
  }
  template<>
  bool BroadCast :: serialize<Pescan::Interface::Escan>( Pescan::Interface::Escan &_p )
  {
    if ( not serialize( _p.filename ) ) return false;
    if ( not serialize( _p.output ) ) return false;
    if ( not serialize( _p.wavefunction_out ) ) return false;
    if ( not serialize( _p.wavefunction_in ) ) return false;
    types::t_unsigned n = _p.method;
    if ( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      switch ( n ) 
      {
        case Pescan::Interface::Escan::NOMET: 
          _p.method = Pescan::Interface::Escan::NOMET; break;
        case Pescan::Interface::Escan::FOLDED_SPECTRUM: 
          _p.method = Pescan::Interface::Escan::FOLDED_SPECTRUM; break;
        default:
        case Pescan::Interface::Escan::ALL_ELECTRON: 
          _p.method = Pescan::Interface::Escan::ALL_ELECTRON; break;
      }
    if ( not serialize( _p.Eref.vbm ) ) return false;
    if ( not serialize( _p.Eref.cbm ) ) return false;
    if ( not serialize( _p.smooth ) ) return false;
    if ( not serialize( _p.kinscal ) ) return false;
    if ( not serialize( _p.nbstates ) ) return false;
    if ( not serialize( _p.niter ) ) return false;
    if ( not serialize( _p.nlines ) ) return false;
    if ( not serialize( _p.tolerance ) ) return false;
    if ( not serialize( &_p.kpoint[0], &_p.kpoint[0]+3 ) ) return false;
    if ( not serialize( _p.scale ) ) return false;
    n = _p.potential;
    if ( not serialize( n ) ) return false;
    if ( stage == COPYING_FROM_HERE )
      switch ( n ) 
      {
        case Pescan::Interface::Escan::NOPOT: 
          _p.potential = Pescan::Interface::Escan::NOPOT; break;
        default:
        case Pescan::Interface::Escan::LOCAL: 
          _p.potential = Pescan::Interface::Escan::LOCAL; break;
        case Pescan::Interface::Escan::NONLOCAL: 
          _p.potential = Pescan::Interface::Escan::NONLOCAL; break;
        case Pescan::Interface::Escan::SPINORBIT: 
          _p.potential = Pescan::Interface::Escan::SPINORBIT; break;
      }
    if ( not serialize( _p.rcut ) ) return false;
    if ( not serialize( _p.launch ) ) return false;
    return serialize_container( _p.spinorbit );
  }
  template<>
  bool BroadCast :: serialize<Pescan::Interface::GenPot>( Pescan::Interface::GenPot &_p )
  {
    if( not serialize( _p.filename ) ) return false;
    if( not serialize( _p.x ) ) return false;
    if( not serialize( _p.y ) ) return false;
    if( not serialize( _p.z ) ) return false;
    if( not serialize( _p.cutoff ) ) return false;
    if( not serialize( _p.output ) ) return false;
    if( not serialize( _p.launch ) ) return false;
    return serialize_container( _p.pseudos );
  }
  template<>
  bool BroadCast :: serialize<Pescan::Interface>( Pescan::Interface &_p )
  {
    if ( not serialize( _p.atom_input ) ) return false;
    if ( not serialize( _p.escan ) ) return false;
    if ( not serialize( _p.genpot ) ) return false;
    types::t_unsigned n = _p.computation;
    if ( not serialize( n ) ) return false;
    if ( stage == COPYING_TO_HERE )
      switch( n )
      {
        case Pescan::Interface::VBM:
          _p.computation = Pescan::Interface::VBM; break;
        default:
        case Pescan::Interface::CBM:
          _p.computation = Pescan::Interface::CBM; break;
      }
    return serialize( _p.dirname );
  }

}
#endif
