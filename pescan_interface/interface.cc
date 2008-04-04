//
//  Version: $Id$
//
#include <fstream>
#include <sstream>
#include <stdexcept>       // std::runtime_error
#ifdef _DIRECTIAGA
#include <unistd.h>
#endif

#include <print/manip.h>
#include <opt/debug.h>

#include <mpi/macros.h>

#include "interface.h"

namespace Pescan
{
  bool Interface :: operator()()
  {
    create_directory();
    create_potential();
    launch_pescan();
    __TRYCODE( return read_result();,
               "Error while reading and/or parsing escan output.\n" ) 
  }

  void Interface :: create_directory()
  {
    std::ostringstream sstr;
    sstr << "mkdir -p " << dirname;
    __DOMPISEQUENTIAL( system( sstr.str().c_str() ); )
    comm->barrier();
  }
  void Interface :: destroy_directory()
  {
    std::ostringstream sstr;
    sstr << "rm -rf " << dirname;
    __DOMPISEQUENTIAL( system( sstr.str().c_str() ); )
  }
  void Interface :: create_potential()
  {
    write_genpot_input();

    std::ostringstream sstr;
    sstr << "cp -u " << atom_input __DIAGA( << "." << comm->rank() ) << " ";
#ifndef _NOLAUNCH
    __IIAGA( sstr <<  genpot.launch << " "; )
#endif
    if( comm->is_root_node() )
    {
      std::vector<std::string> :: const_iterator i_str = genpot.pseudos.begin();
      std::vector<std::string> :: const_iterator i_str_end = genpot.pseudos.end();
      for(; i_str != i_str_end; ++i_str)  sstr << *i_str << " ";
    }
    sstr << " " << dirname;
    // Makes sure that procs don't simultaneously access the same file.
    __DOMPISEQUENTIAL( system( sstr.str().c_str() ); )

    
#ifndef _NOLAUNCH
    chdir( dirname.c_str() );
    __IIAGA( 
      sstr.str("");
      sstr << "./" << Print::StripDir(genpot.launch);
      system(sstr.str().c_str()); 
    )
    __DIAGA( 
      int __rank = comm->rank();
      MPI_Comm __commC = (MPI_Comm) *( (MPI::Comm*) comm->get() ) ;
      MPI_Fint __commF = MPI_Comm_c2f( __commC );
      FC_FUNC_(iaga_call_genpot, IAGA_CALL_GENPOT)( &__commF, &__rank );
    )
    chdir( ".." );
#endif
  }
  types::t_real Interface :: launch_pescan()
  {
    std::ostringstream sstr;
    write_escan_input();

    sstr << "cp -u " << maskr << " "; 
#ifndef _NOLAUNCH
    __IIAGA(sstr << escan.launch << " ";)
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
    __DOMPISEQUENTIAL( system( sstr.str().c_str() ); ) 

#ifndef _NOLAUNCH
    chdir( dirname.c_str() );
    __IIAGA(
      std::string output = Print::StripEdges(escan.output);
      sstr.str("");
      sstr << " ./" << Print::StripDir(escan.launch) << " > " << output;
      system(sstr.str().c_str());
    )
    __DIAGA( FC_FUNC_(iaga_call_escan, IAGA_CALL_ESCAN)( &escan.nbstates ); )
    chdir( ".." );
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
    __DOASSERT( not child, "No <GenPot> tag found on input.\nAborting\n" )
    __DOASSERT( not child->Attribute("x"), "x attributes are missing in <GenPot>.\n")
    __DOASSERT( not child->Attribute("y"), "y attributes are missing in <GenPot>.\n")
    __DOASSERT( not child->Attribute("z"), "z attributes are missing in <GenPot>.\n")
    __DOASSERT( not child->Attribute("cutoff"),
                "cutoff attributes are missing in <GenPot>.\n")
    __DOASSERT( not child->FirstChildElement("Pseudo"),
                "Could not find <Pseudo/> in <GenPot>\nAborting\n")

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
    __DOASSERT( not genpot.check(),
                "Insufficient or Incorrect Genpot input.\nAborting\n" )


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
    __DOASSERT( not child, "Could not find <Hamiltonian> tag on input.\n")
    __DOASSERT( not child->Attribute("kinscal"),
                "kinscal attribute is missing in <Hamiltonian> tag.\n")
    __DOASSERT( not child->Attribute("realcutoff"),
                "realcutoff attribute is missing in <Hamiltonian> tag.\n")
    __DOASSERT( not child->Attribute("potential"),
                "potential attribute is missing in <Hamiltonian> tag.\n")
    if ( child->Attribute("launch") )
      escan.launch = child->Attribute("launch");

    __TRYCODE( check_existence();, "Some files and/or programs are missing.\n" )

    child->Attribute("kinscal", &escan.kinscal);
    if( child->Attribute("smooth") )
      child->Attribute("smooth", &escan.smooth);
    if( child->Attribute("potential") )
    {
      child->Attribute("potential", &j);
      switch( j )
      {
        default:
        case Escan::NOPOT: __THROW_ERROR( "Error, incorrect escan potential\n" ) break;
        case Escan::LOCAL: escan.potential = Escan::LOCAL; break;
        case Escan::NONLOCAL: escan.potential = Escan::NONLOCAL; break;
        case Escan::SPINORBIT: escan.potential = Escan::SPINORBIT; break;
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
    if( not child ) return true;

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
    std::ostringstream sstr;
    sstr << Print::StripEdges(dirname)
         << "/" << Print::StripEdges(genpot.filename)
         __DIAGA( << "." << comm->rank() );
    file.open( sstr.str().c_str(), std::ios_base::out|std::ios_base::trunc ); 

    __DOASSERT( file.bad() or ( not file.is_open() ),
                   "Could not open file " << sstr.str()
                << " for writing.\nAborting.\n" )

    file << Print::StripDir(atom_input) << "\n"
         << genpot.x << " " 
         << genpot.y << " " 
         << genpot.z << "\n"
         << genpot.x << " " 
         << genpot.y << " " 
         << genpot.z << "\n" << " 0 0 0\n" 
         << genpot.cutoff << "\n"
         << genpot.pseudos.size() << std::endl;
    std::vector< std::string > :: const_iterator i_str = genpot.pseudos.begin();
    std::vector< std::string > :: const_iterator i_str_end = genpot.pseudos.end();
    for(; i_str != i_str_end; ++i_str )
     file << Print::StripDir(*i_str) << "\n";
    file.flush();
    file.close();
  }

  void Interface::write_escan_input()
  {
    std::ofstream file;
    std::ostringstream sstr;
    sstr << Print::StripEdges(dirname) << "/"
         __IIAGA( << Print::StripEdges(escan.filename) )
         __DIAGA( << "escan_input." << comm->rank() );
    std::string name = sstr.str();
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 

    __DOASSERT( file.bad() or ( not file.is_open() ),
                   "Could not open file " << (std::string) name
                << " for writing.\nAborting.\n" )

    file << "1 " << Print::StripDir(dirname, genpot.output) 
                 __DIAGA( << "." << comm->rank() ) << "\n"
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
      atat::rVector3d k = escan.kpoint;
      file << "11 1 " << k << " " 
                      << std::setw(12) << std::setprecision(8) 
                      << escan.scale / Physics::a0("A") <<  "\n";
    }
    file << "12 " << escan.potential << "\n"
         << "13 " << atom_input __DIAGA( << "." << comm->rank() ) << "\n"
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
#ifndef _DIRECTIAGA
    std::ifstream file;
    std::ostringstream sstr;
    sstr << dirname; 
    sstr << "/" << Print::StripEdges(escan.output);
    std::string name = sstr.str();
    file.open( name.c_str(), std::ios_base::in ); 
    __DOASSERT( not file.is_open(), 
                "Could not open file " << (std::string) name << "\n" )

    char cline[256];
    std::string line("");
    while (     ( not file.eof() )
            and line.find("FINAL eigen") == std::string::npos ) 
    {
      file.getline( cline, 256 );
      line = cline;
    }
    __DOASSERT( file.eof(),
                   "Reached end of " << name
                << " without encoutering eigenvalues\n" )

    eigenvalues.clear(); eigenvalues.reserve(escan.nbstates);
    types::t_unsigned u(0);
    types::t_real eig;
    for(; u < escan.nbstates; ++u )
    {
      if(  file.bad() ) break;
      file >> eig;
      eigenvalues.push_back( eig );
    }
#endif

    __DIAGA(
      double values[ escan.nbstates ];
      eigenvalues.resize( escan.nbstates );
      FC_FUNC_(iaga_get_eigenvalues, IAGA_GET_EIGENVALUES)( values, &escan.nbstates );
      comm->get()->Bcast( values, escan.nbstates, MPI::DOUBLE, 0 );
      std::copy( values, values + escan.nbstates, eigenvalues.begin() );
    )

    __IIAGA( __DOASSERT( u != escan.nbstates,
                            "Found " << u << " eigenvalues in " << name
                         << " where " << escan.nbstates
                         << " were expected.\n" ) )
    return true;
  }
               
}
