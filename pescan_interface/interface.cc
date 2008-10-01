//
//  Version: $Id$
//
#include <fstream>
#include <sstream>
#include <stdexcept>       // std::runtime_error
#ifdef _DIRECTIAGA
# include <unistd.h>
#endif
#include <boost/mpi/collectives.hpp>

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
    namespace bfs = boost::filesystem;
    const t_Path newpath(  InitialPath::path() / dirname );
    __DIAGA
    ( 
      for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i )
        if( MPI_COMM.rank() == i )
    )
    if( not bfs::exists( newpath ) ) bfs::create_directory( newpath );
  }
  void Interface :: destroy_directory()
  {
    namespace bfs = boost::filesystem;
    const t_Path newpath(  InitialPath::path() / dirname );
    __DIAGA
    ( 
     for( types::t_int i( 1 ); i < MPI_COMM.size(); ++i )
       if( MPI_COMM.rank() == i )
    )
    if( bfs::exists( newpath ) ) bfs::remove_directory( newpath );
  }
  void Interface :: create_potential()
  {
    namespace bfs = boost::filesystem;
    write_genpot_input();
    
    const t_Path( atom_input __DIAGASUFFIX )
    bfs::copy_file( InitialPath::path()/orig, InitialPath::path()/dirname );

#   ifndef _NOLAUNCH
      __IIAGA( bfs::copy_file( genpot.launch, InitialPath::path()/dirname ); )
#   endif
    __DIAGA
    (
      for( types::t_int i( 1 ); i < MPI_COMM.size(); ++i )
        if( MPI_COMM.rank() == i )
    )
    {
      std::vector<t_Path> :: const_iterator i_str = genpot.pseudos.begin();
      std::vector<t_Path> :: const_iterator i_str_end = genpot.pseudos.end();
      for(; i_str != i_str_end; ++i_str) 
        bfs::copy_file( *i_str, InitialPath::path()/dirname );
    }

    
#   ifndef _NOLAUNCH
      chdir( ( InitialPath::path() / dirname ).c_str() );
      __IIAGA(  system( "./" + genpot.launch.filename().string() ) );
      __DIAGA
      ( 
        int __rank = MPI_COMM.rank();
        MPI_Comm __commC = (MPI_Comm) ( MPI_COMM ) ;
        MPI_Fint __commF = MPI_Comm_c2f( __commC );
        FC_FUNC_(iaga_call_genpot, IAGA_CALL_GENPOT)( &__commF, &__rank );
      )
      chdir( InitialPath::path().c_str() );
#   endif
  }
  types::t_real Interface :: launch_pescan()
  {
    write_escan_input();

    __DIAGA
    ( 
      for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i )
        if( MPI_COMM.rank() == i )
    )
    bfs::copy_file( maskr, InitialPath::path()/dirname );
#   ifndef _NOLAUNCH
      __IIAGA( bfs::copy_file( escan.launch, InitialPath::path()/dirname ); )
#   endif
    std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
    std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
    for( ; i_so != i_so_end; ++i_so)
      __DIAGA
      ( 
        for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i )
          if( MPI_COMM.rank() == i )
      )
      bfs::copy_file( i_so->filename, InitialPath::path()/dirname );
    
#   ifndef _NOLAUNCH
      chdir( ( InitialPath::path() / dirname ).c_str() );
      __IIAGA(  system( "./" + escan.launch.filename().string() ) );
      __DIAGA( FC_FUNC_(iaga_call_escan, IAGA_CALL_ESCAN)( &escan.nbstates ); )
      chdir( InitialPath::path().c_str() );
#   endif
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
      maskr = child->Attribute("filename");

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
        genpot.pseudos.push_back( child->Attribute("filename") );
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
    const t_Path orig(   opt::InitialPath::path() 
                       / dirname / genpot.filename __DIAGASUFFIX );
    file.open( orig.string(), std::ios_base::out|std::ios_base::trunc ); 

    __DOASSERT( file.bad() or ( not file.is_open() ),
                   "Could not open file " << sstr.str()
                << " for writing.\nAborting.\n" )

    file << atom_input << "\n"
         << genpot.x << " " 
         << genpot.y << " " 
         << genpot.z << "\n"
         << genpot.x << " " 
         << genpot.y << " " 
         << genpot.z << "\n" << " 0 0 0\n" 
         << genpot.cutoff << "\n"
         << genpot.pseudos.size() << std::endl;
    std::vector< t_Path > :: const_iterator i_str = genpot.pseudos.begin();
    std::vector< t_Path > :: const_iterator i_str_end = genpot.pseudos.end();
    for(; i_str != i_str_end; ++i_str )
     file << Print::StripDir(*i_str) << "\n";
    file.flush();
    file.close();
  }

  void Interface::write_escan_input()
  {
    std::ofstream file;
    const t_Path orig
    ( 
        opt::InitialPath::path() / dirname 
      / __IIAGA( genpot.filename )
        __DIAGA( "escan_input." ) __DIAGASUFFIX
    );
    std::ostringstream sstr;
    std::string name = orig.string();
    file.open( name.c_str(), std::ios_base::out|std::ios_base::trunc ); 

    __DOASSERT( file.bad() or ( not file.is_open() ),
                   "Could not open file " << (std::string) name
                << " for writing.\nAborting.\n" )

    file << "1 " << ( genpot.output __DIAGASUFFIX ) 
         << "2 " << escan.wavefunction_out << "\n"
         << "3 " << escan.method << "\n"
         << "4 " << escan.Eref << " " << genpot.cutoff << " "
                 << escan.smooth << " " << escan.kinscal << "\n"
         << "5 " << escan.nbstates << "\n"
         << "6 " << escan.niter << " " << escan.nlines
                 << " " << escan.tolerance << "\n";
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
         << "13 " << ( atom_input __DIAGASUFFIX ) << "\n"
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
#   ifndef _NOLAUNCH
#     ifndef _DIRECTIAGA
        std::ifstream file;
        std::ostringstream sstr;
        const t_Path orig( opt::InitialPath::path() / dirname / escan.output );
        std::string name = orig.string();
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
#     endif
     
      __DIAGA(
        double values[ escan.nbstates ];
        eigenvalues.resize( escan.nbstates );
        FC_FUNC_(iaga_get_eigenvalues, IAGA_GET_EIGENVALUES)(values,&escan.nbstates);
        boost::mpi::broadcast( MPI_COMM, values, escan.nbstates, 0 );
        std::copy( values, values + escan.nbstates, eigenvalues.begin() );
          Print::out << "Eigenvalues: " << Print::endl;
        for( types :: t_int i = 0; i < escan.nbstates; i++) 
          Print::out << eigenvalues[i] << " ---> " << i << Print::endl;
      )
     
      __IIAGA( __DOASSERT( u != escan.nbstates,
                              "Found " << u << " eigenvalues in " << name
                           << " where " << escan.nbstates
                           << " were expected.\n" ) )
#   endif
    return true;
  }
               
}
