//
//  Version: $Id$
//
#include <fstream>
#include <sstream>
#include <stdexcept>       // std::runtime_error

#include <print/manip.h>

#include "interface.h"

namespace Pescan
{
  bool Interface :: operator()()
  {
    create_directory();
    create_potential();
    launch_pescan();
    return read_result();
  }

  void Interface :: create_directory()
  {
    std::ostringstream sstr;
    sstr << "mkdir " << dirname;
    system( sstr.str().c_str() );
  }
  void Interface :: destroy_directory()
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
  types::t_real Interface :: launch_pescan()
  {
    std::ostringstream sstr;
    write_escan_input();

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

    // This whole section tries to find a <Functional type="escan"> tag
    // in _element or its child
  const TiXmlElement* Interface :: find_node (const TiXmlElement &_node )
  {
    const TiXmlElement *parent;
    std::string str;

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
      if ( str.compare("escan" ) == 0 ) break;
      parent = parent->NextSiblingElement("Functional");
    }

    return parent;
  }

  bool Interface :: Load_ (const TiXmlElement &_node )
  {
    const TiXmlElement *child;

    // lester needs to load "module fftw" when in interactive multi-proc mode
    if ( _node.Attribute("system") )
    {
      std::string str = _node.Attribute("system");
      system( str.c_str() );
    }

    do_destroy_dir = true;
    if( _node.Attribute("keepdirs") ) do_destroy_dir = false;

    child = _node.FirstChildElement("GenPot");
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


    child = _node.FirstChildElement("Wavefunctions");
    if( child and child->Attribute("in") )
      escan.wavefunction_in = child->Attribute("in");
    if( child and child->Attribute("out") )
      escan.wavefunction_out = child->Attribute("out");
    types::t_int j;
    if( _node.Attribute("method", &j) )
      escan.method = ( j == 1 ) ? FOLDED_SPECTRUM: ALL_ELECTRON;
    child = _node.FirstChildElement("Reference");
    if( child and child->Attribute("value") )
      child->Attribute("value", &escan.Eref);

    child = _node.FirstChildElement("Hamiltonian");
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
    child = _node.FirstChildElement("Minimizer");
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

  void Interface::write_escan_input() 
  {
    types::t_unsigned nbwfn = escan.nbstates; 
    std::ofstream file;
    std::string name = Print::StripEdges(dirname) + "/" + Print::StripEdges(escan.filename);
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    file << "1 " << Print::StripDir(dirname, genpot.output) << "\n"
         << "2 " << escan.wavefunction_out << "\n"
         << "3 " << escan.method << "\n"
         << "4 " << escan.Eref << " " << genpot.cutoff << " "
                 << escan.smooth << " " << escan.kinscal << "\n"
         << "5 " << nbwfn << "\n"
         << "6 " << escan.niter << " " << escan.nlines << " " << escan.tolerance << "\n";
    if( do_input_wavefunctions )
    {
      file << "7 " << escan.nbstates
           << "\n8 ";
      for( types::t_unsigned u=1; u < escan.nbstates; ++u )
        file << u << ", ";
      file << escan.nbstates 
           << "\n9 " << escan.wavefunction_in << "\n";
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

  bool Interface :: read_result()
  {
#ifdef _NOLAUNCH
    return true;
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

    eigenvalues.clear(); eigenvalues.reserve(escan.nbstates);
    types::t_unsigned u(0);
    for(; u < escan.nbstates; ++u )
    {
      types::t_real eig;
      if( not file.eof() ) break;
      if( (file >> eig).good() ) break;
      eigenvalues.push_back( eig );
    }
    if( u == escan.nbstates ) return true;
    
    std::cerr << "Found " << u << " eigenvalues when " 
              << escan.nbstates << " expected " << std::endl;
    return false;
  }
               

}

#ifdef _MPI

namespace mpi
{
  template<>
  bool BroadCast :: serialize<Pescan::Interface::SpinOrbit>( Pescan::Interface::SpinOrbit &_sp )
  {
    return     serialize( _sp.filename ) 
           and serialize( _sp.izz )
           and serialize( _sp.s ) 
           and serialize( _sp.p ) 
           and serialize( _sp.d ) 
           and serialize( _sp.pnl ) 
           and serialize( _sp.dnl );
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
        case Pescan::Interface::NOMET: 
          _p.method = Pescan::Interface::NOMET; break;
        case Pescan::Interface::FOLDED_SPECTRUM: 
          _p.method = Pescan::Interface::FOLDED_SPECTRUM; break;
        default:
        case Pescan::Interface::ALL_ELECTRON: 
          _p.method = Pescan::Interface::ALL_ELECTRON; break;
      }
    if ( not serialize( _p.Eref.vbm ) ) return false;
    if ( not serialize( _p.Eref.cbm ) ) return false;
    if ( not serialize( _p.smooth ) ) return false;
    if ( not serialize( _p.kinscal ) ) return false;
    if ( not serialize( _p.escan.nbstates ) ) return false;
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
    return     serialize( _p.filename ) 
           and serialize( _p.x ) 
           and serialize( _p.y ) 
           and serialize( _p.z ) 
           and serialize( _p.cutoff ) 
           and serialize( _p.output ) 
           and serialize( _p.launch ) 
           and serialize_container( _p.pseudos );
  }
  template<>
  bool BroadCast :: serialize<Pescan::Interface>( Pescan::Interface &_p )
  {
    return     serialize( _p.atom_input )
           and serialize( _p.escan )       
           and serialize( _p.genpot )    
           and serialize( _p.dirname )
           and serialize_container(_p.eigenvalues );
  }

}
#endif
