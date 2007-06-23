#include "interface.h"
#include<mpi/mpi_object.h>
#include <sstream>

#include<physics/physics.h>

namespace Pescan
{
  types::t_real Interface :: operator()( Ising_CE::Structure &_str)
  {
    if ( escan.method == Escan::FOLDED_SPECTRUM )
      computation = CBM;
      
    escan.scale = _str.scale;
    create_directory();
    create_potential();
    launch_pescan( _str );
    types::t_real result = read_result( _str );

    if ( escan.method == Escan::ALL_ELECTRON )
    {
      destroy_directory();
      return result;
    }
    band_edge.second = result;

    computation = VBM;
    launch_pescan( _str );
    band_edge.first = read_result( _str );
    destroy_directory();
    return result - band_edge.first;
  }
  void Interface :: create_directory()
  {
    std::ostringstream sstr;
    sstr << "escan.";
#ifdef _MPI
    sstr << mpi::main.rank();
#endif
    dirname = sstr.str();

    sstr.str("");
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
    std::ostringstream sstr;
    sstr << "cp " << atom_input << " " << dirname;
    system( sstr.str().c_str() );

    write_genpot_input();
    sstr.str(""); 
    sstr << "mv " << genpot.filename << " " << dirname;
    system( sstr.str().c_str() );

    sstr.str("");
    sstr << "cp " << genpot.launch << " ";
    std::vector<std::string> :: const_iterator i_str = genpot.pseudos.begin();
    std::vector<std::string> :: const_iterator i_str_end = genpot.pseudos.end();
    for(; i_str != i_str_end; ++i_str)
      sstr << *i_str << " ";
    sstr << " " << dirname;
    system( sstr.str().c_str() );

    
    sstr.str("");
    sstr << "cd " << dirname << ";" << genpot.launch;
    system(sstr.str().c_str());
  }
  types::t_real Interface :: launch_pescan( Ising_CE::Structure &_str )
  {
    std::ostringstream sstr;
    write_escan_input( _str );
    sstr << "mv " << escan.filename << " " << dirname;
    system( sstr.str().c_str() );

    sstr.str("");
    sstr << "cp " << escan.launch << " ";
    sstr << " maskr ";
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

    sstr.str("");
    sstr << "cd " << dirname << ";" << escan.launch << " > " << escan.output;
    system(sstr.str().c_str());
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
        child->Attribute("VBM", &escan.Eref.first);
      if( child->Attribute("CBM") )
        child->Attribute("CBM", &escan.Eref.second);
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

  }

  void Interface::write_genpot_input()
  {
    std::ofstream file;
    file.open( genpot.filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 

    file << atom_input << std::endl 
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
     file << *i_str << std::endl;
    file.flush();
    file.close();
  }

  void Interface::write_escan_input( Ising_CE::Structure &_str ) 
  {
    std::ofstream file;
    file.open( escan.filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    file << "1 " << genpot.output << std::endl
         << "2 " << escan.wavefunction << std::endl
         << "3 " << escan.method << std::endl;
    switch (computation)
    { 
      case VBM: file << "4 " << escan.Eref.first << " "; break;
      case CBM: file << "4 " << escan.Eref.second << " "; break;
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
      file << "5 " << n + escan.nbstates << std::endl;
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
    sstr << dirname << "/" << escan.output;
    file.open( sstr.str().c_str(), std::ios_base::in ); 
    char cline[256];
    std::string line("");
    while (     ( not file.eof() )
            and line.find("FINAL eigen") == std::string::npos ) 
    {
      file.getline( cline, 256 );
      line = cline;
    }

    std::vector<types::t_real> eigenvalues;
    types :: t_real eigenvalue;
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
      band_edge.first = eigenvalues[n-1];
      band_edge.second = eigenvalues[n];
      return eigenvalues[n] - eigenvalues[n-1];
    }

    types :: t_real reference = 0;
    if ( escan.method == Escan::FOLDED_SPECTRUM )
      switch (computation)
      { 
        case VBM: reference = escan.Eref.first; break;
        case CBM: reference = escan.Eref.second; break;
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
    if ( not serialize( _p.Eref.first ) ) return false;
    if ( not serialize( _p.Eref.second ) ) return false;
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
