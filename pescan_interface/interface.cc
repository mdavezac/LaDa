//
//  Version: $Id$
//
#include "interface.h"
#include<mpi/mpi_object.h>
#include <sstream>
#include <stdlib.h>

#include<physics/physics.h>

std::string StripDir( std::string _string )
{
  _string = StripEdges( _string );
  size_t t = _string.find_last_of("/");
  if ( t == std::string::npos )
    return _string;
  
  return _string.erase( 0, t + 1 );
}
std::string StripDir( const std::string &_dir, const std::string &_str )
{
  std::string dir = reformat_home(_dir);
  std::string str = reformat_home(_str);
  if( str.find( dir ) == std::string::npos )
    return str;
  
  str.erase(0, dir.length() );
  if ( _str.find_first_of(" /\t\n") )
    str.erase(0, str.find_first_not_of(" /\t\n"));
  return str; 
}
std::string StripEdges( std::string _string )
{
  if ( _string.find_first_of(" \t\n") )
    _string.erase(0, _string.find_first_not_of(" \t\n") );
  size_t l = _string.length();
  size_t t = _string.find_last_not_of(" \t\n") + 1;
  if ( t > l ) return _string;
  _string.erase(t, l );
  return _string;
}

std::string reformat_home ( std::string _str )
{
  // Trim leading and tailing spaces
  _str = StripEdges(_str);
  if ( _str[0] != '~' ) return _str;
  std::string home="";
  if ( getenv("HOME") ) home = getenv("HOME");
  else if ( getenv("home") ) home = getenv("home");
  else return _str;
  if ( _str.find_first_of("/") )
    _str.erase(0, _str.find_first_of("/") );
  _str.insert(0, home );
  return _str;
}

namespace Pescan
{
  types::t_real Interface :: operator()( Ising_CE::Structure &_str)
  {
#ifdef _NOLAUNCH
    bands.vbm = random() * 5;
    bands.cbm = bands.vbm + random() * 5;
    return bands.gap();
#endif
    if ( escan.method == Escan::FOLDED_SPECTRUM )
      computation = CBM;
      
    std::string olddirname = dirname;
    escan.scale = _str.scale;
    if ( escan.method == Escan::FOLDED_SPECTRUM and not do_destroy_dir )
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
      destroy_directory();
      escan.Eref = bands;
      return result;
    }
    bands.cbm = result;

    computation = VBM;
    destroy_directory();
    if ( escan.method == Escan::FOLDED_SPECTRUM and not do_destroy_dir )
      dirname = olddirname + ( computation == CBM ?  "/cbm": "/vbm" );
    create_directory();
    create_potential();
    launch_pescan( _str );
    bands.vbm = read_result( _str );
    destroy_directory();
    return bands.gap();
  }
  void Interface :: create_directory()
  {
    if ( do_destroy_dir ) 
    {
      std::ostringstream sstr;
      sstr << "escan";
#ifdef _MPI
      sstr << "." << mpi::main.rank();
#endif
      dirname = sstr.str();
    }

    std::ostringstream sstr;
    sstr << "mkdir " << dirname;
    system( sstr.str().c_str() );
  }
  void Interface :: destroy_directory()
  {
    if ( not do_destroy_dir ) return;
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
    sstr << "cd " << dirname << "; ./" << StripDir(genpot.launch);
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

    std::string output = StripEdges(escan.output);
    std::cout << "Method " << (escan.method == Escan::FOLDED_SPECTRUM ? "Folded": "All Electron") << std::endl;
    if ( escan.method == Escan::FOLDED_SPECTRUM )
      output += ( computation == VBM ) ? ".vbm": ".cbm";
    sstr.str("");
    sstr << "cd " << dirname << "; ./" << StripDir(escan.launch) << " > " << output;
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


    child = parent->FirstChildElement("wavefunction");
    if( child and child->Attribute("name") )
      escan.wavefunction = child->Attribute("name");
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
    std::string name = StripEdges(dirname) + "/" + StripEdges(genpot.filename);
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 

    file << StripDir(atom_input) << std::endl 
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
     file << StripDir(*i_str) << std::endl;
    file.flush();
    file.close();
  }

  void Interface::write_escan_input( Ising_CE::Structure &_str ) 
  {
    std::ofstream file;
    std::string name = StripEdges(dirname) + "/" + StripEdges(escan.filename);
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    file << "1 " << StripDir(dirname, genpot.output) << std::endl
         << "2 " << escan.wavefunction << std::endl
         << "3 " << escan.method << std::endl;
    switch (computation)
    { 
      case VBM: file << "4 " << escan.Eref.vbm << " "; break;
      case CBM: file << "4 " << escan.Eref.cbm << " "; break;
    }
    file << genpot.cutoff << " " 
         << escan.smooth << " "
         << escan.kinscal << std::endl;

    if( escan.method == Escan::FOLDED_SPECTRUM )
      file << "5 " << escan.nbstates << std::endl;
    else
    {
      Ising_CE::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
      Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
      types::t_unsigned n = 0; 
      for(; i_atom != i_atom_end; ++i_atom)
      {
        Ising_CE::StrAtom atom; 
        _str.lattice->convert_Atom_to_StrAtom( *i_atom, atom );
        n += physics::atoms::Charge( atom.type );
      }
      if (    escan.potential != Escan::SPINORBIT
           or atat::norm2(escan.kpoint) < types::tolerance )
        n /= 2;
      file << "5 " << n + 3 << std::endl;
    }
    file << "6 " << escan.niter << " " << escan.nlines << " " << escan.tolerance << std::endl
         << "7 0" << std::endl << "8 0" << std::endl << "9 wg.cbm.in" << std::endl
         << "10 0 4 4 4 graph.R" << std::endl;
    if ( atat::norm2( escan.kpoint ) < types::tolerance )
     file << "11 0 0 0 0 ";
    else 
     file << "11 1 " << escan.kpoint[0] << " " << escan.kpoint[1] << " " << escan.kpoint[2] << " ";
    file << escan.scale << std::endl
         << "12 " << escan.potential << std::endl
         << "13 " << atom_input << std::endl
         << "14 " << escan.rcut << std::endl
         << "15 " << escan.spinorbit.size() << std::endl;
    std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
    std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
    for(types::t_unsigned i=16; i_so != i_so_end; ++i_so, ++i )
    {
      file << i << " " << i_so->filename << " " 
           << i_so->izz << " " 
           << i_so->s << " " << i_so->p << " " << i_so->d << " " 
           << i_so->pnl << " " << i_so->dnl << std::endl;
    }
    file.flush();
    file.close();
  }

  types::t_real Interface :: read_result( Ising_CE::Structure &_str )
  {
    std::ifstream file;
    std::ostringstream sstr;
    sstr << dirname << "/" << StripEdges(escan.output);
    if ( escan.method == Escan::FOLDED_SPECTRUM )
      sstr << ( ( computation == VBM ) ? ".vbm": ".cbm" );
    file.open( sstr.str().c_str(), std::ios_base::in ); 
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
        n += physics::atoms::Charge( atom.type );
      }
      if (    escan.potential != Escan::SPINORBIT
           or atat::norm2(escan.kpoint) < types::tolerance )
        n /= 2;
      // watch out! C arrays start at index = 0
      if ( eigenvalues.size() < n )
      {
        std::cerr << "Error, not enough states were computed to determine band gap " << n << std::endl;
        return -1.0;
      }
      bands.vbm = eigenvalues[n-1];
      bands.cbm = eigenvalues[n];
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
    if ( not serialize( _p.wavefunction ) ) return false;
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
