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

    sstr << "cp " << maskr << " "; 
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

    child = _node.FirstChildElement("Maskr");
    if( child and child->Attribute("filename") )
      maskr = Print::reformat_home(child->Attribute("filename"));

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
        genpot.pseudos.push_back( Print::reformat_home(child->Attribute("filename")) );
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
        so.filename = Print::reformat_home(child->Attribute("filename"));
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
          child->Attribute("dnl", &so.dnl);
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
    std::ofstream file;
    std::string name = Print::StripEdges(dirname) + "/" + Print::StripEdges(escan.filename);
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    if( file.bad() or ( not file.is_open() ) )
    {
      std::ostringstream error;
      error << __FILE__ << ", line: " << __LINE__ << ":\n"
            <<  " Could not open file " << name << " for writing.\n"
            <<  " Aborting." << std::endl;
      throw std::runtime_error( error.str() );
    }
    file << "1 " << Print::StripDir(dirname, genpot.output) << "\n"
         << "2 " << escan.wavefunction_out << "\n"
         << "3 " << escan.method << "\n"
         << "4 " << escan.Eref << " " << genpot.cutoff << " "
                 << escan.smooth << " " << escan.kinscal << "\n"
         << "5 " << escan.nbstates << "\n"
         << "6 " << escan.niter << " " << escan.nlines << " " << escan.tolerance << "\n";
    if( escan.read_in.size() )
    {
      file << "7 " << escan.read_in.size() << "\n"
           << "8 ";
      std::vector<types::t_unsigned> :: const_iterator i_r = escan.read_in.begin();
      std::vector<types::t_unsigned> :: const_iterator i_r_end = escan.read_in.end();
      file << *i_r;
      for(--i_r_end; i_r != i_r_end; ++i_r)
        file << ", " << *i_r;
      file << "\n" << "9 " << escan.wavefunction_in << "\n"
           << "10 0 4 4 4 graph.R\n";
    }
    else  file << "7 0\n8 0\n9 dummy\n10 0 4 4 4 graph.R\n";

    if ( atat::norm2( escan.kpoint ) < types::tolerance ) file << "11 0 0 0 0 0\n";
    else
    {
      atat::rVector3d k = escan.kpoint * 2.0 * Math::pi * Physics::a0("A") / escan.scale;
      file << "11 1 " << k << " " 
                      << std::setw(12) << std::setprecision(8) 
                      << escan.scale / Physics::a0("A") <<  "\n";
    }
    file << "12 " << escan.potential << "\n"
         << "13 " << atom_input << "\n"
         << "14 " << escan.rcut << "\n"
         << "15 " << escan.spinorbit.size() << "\n";
    std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
    std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
    for(types::t_unsigned i=16; i_so != i_so_end; ++i_so, ++i )
      file << i << " " << Print::StripDir(i_so->filename) << " " 
           << i_so->izz << " " 
           << i_so->s << " " << i_so->p << " " << i_so->d << " " 
           << i_so->pnl << " " << i_so->dnl << "\n";
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
      sstr << __FILE__ << ", line " << __LINE__ << ":\n"
           << "Could not open file "
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
    types::t_real eig;
    for(; u < escan.nbstates; ++u )
    {
      if(  file.bad() ) break;
      file >> eig;
      eigenvalues.push_back( eig );
    }
    if( u == escan.nbstates ) return true;
    
    std::cerr << __FILE__ << ", line: " << __LINE__ << "\n"
              << "Found " << u << " eigenvalues when " 
              << escan.nbstates << " expected " << std::endl;
    return false;
  }
               
#ifdef _MPI 
  bool Interface::SpinOrbit::broadcast( mpi::BroadCast &_bc )
  {
    return     _bc.serialize( filename ) 
           and _bc.serialize( izz )
           and _bc.serialize( s ) 
           and _bc.serialize( p ) 
           and _bc.serialize( d ) 
           and _bc.serialize( pnl ) 
           and _bc.serialize( dnl );
  }
  bool Interface::Escan::broadcast( mpi::BroadCast &_bc )
  {
    if ( not _bc.serialize( filename ) ) return false;
    if ( not _bc.serialize( output ) ) return false;
    if ( not _bc.serialize( wavefunction_out ) ) return false;
    if ( not _bc.serialize( wavefunction_in ) ) return false;
    if ( not _bc.serialize_container( read_in ) ) return false;
    types::t_unsigned n = method;
    if ( not _bc.serialize( n ) ) return false;
    if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
      switch ( n ) 
      {
        case Interface::NOMET: 
          method = Interface::NOMET; break;
        case Interface::FOLDED_SPECTRUM: 
          method = Interface::FOLDED_SPECTRUM; break;
        default:
        case Interface::ALL_ELECTRON: 
          method = Interface::ALL_ELECTRON; break;
      }
    if ( not _bc.serialize( Eref ) ) return false;
    if ( not _bc.serialize( smooth ) ) return false;
    if ( not _bc.serialize( kinscal ) ) return false;
    if ( not _bc.serialize( nbstates ) ) return false;
    if ( not _bc.serialize( niter ) ) return false;
    if ( not _bc.serialize( nlines ) ) return false;
    if ( not _bc.serialize( tolerance ) ) return false;
    if ( not _bc.serialize( &kpoint[0], &kpoint[0]+3 ) ) return false;
    if ( not _bc.serialize( scale ) ) return false;
    n = potential;
    if ( not _bc.serialize( n ) ) return false;
    if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
      switch ( n ) 
      {
        case Interface::Escan::NOPOT: 
          potential = Interface::Escan::NOPOT; break;
        default:
        case Interface::Escan::LOCAL: 
          potential = Interface::Escan::LOCAL; break;
        case Interface::Escan::NONLOCAL: 
          potential = Interface::Escan::NONLOCAL; break;
        case Interface::Escan::SPINORBIT: 
          potential = Interface::Escan::SPINORBIT; break;
      }
    if ( not _bc.serialize( rcut ) ) return false;
    if ( not _bc.serialize( launch ) ) return false;
    n = spinorbit.size();
    if ( not _bc.serialize(n) ) return false;
    if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
      spinorbit.resize(n);
    std::vector< SpinOrbit > :: iterator i_so = spinorbit.begin();
    std::vector< SpinOrbit > :: iterator i_so_end = spinorbit.end();
    for(; i_so != i_so_end; ++i_so )
      if( not i_so->broadcast( _bc ) ) return false;

    return true;
  }
  bool Interface::GenPot::broadcast( mpi::BroadCast& _bc )
  {
    return     _bc.serialize( filename ) 
           and _bc.serialize( x ) 
           and _bc.serialize( y ) 
           and _bc.serialize( z ) 
           and _bc.serialize( cutoff ) 
           and _bc.serialize( output ) 
           and _bc.serialize( launch ) 
           and _bc.serialize_container( pseudos );
  }
#endif // _MPI 

}

#ifdef _MPI

namespace mpi
{
  template<>
  bool BroadCast :: serialize<Pescan::Interface>( Pescan::Interface &_p )
  {
    return     serialize( _p.atom_input )
           and _p.escan.broadcast( *this )       
           and _p.genpot.broadcast( *this )    
           and serialize( _p.dirname )
           and serialize_container(_p.eigenvalues );
  }

}
#endif
