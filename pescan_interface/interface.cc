//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sstream>
#include <fstream>
#include <stdexcept>       // std::runtime_error
#ifdef _DIRECTIAGA
# include <unistd.h>
# include <boost/mpi/collectives.hpp>
#endif
#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/lexical_cast.hpp>

#include <opt/debug.h>
#include <opt/initial_path.h>
#include <mpi/macros.h>
#include <opt/initial_path.h>

#include "interface.h"


#ifdef _DIRECTIAGA
  // declares fortran interface
  //! \cond
  extern "C"
  {
    void FC_FUNC_(getvlarg, GETVLARG)();
    void FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( MPI_Fint * );
    void FC_FUNC_(iaga_call_escan, IAGA_CALL_ESCAN)( int* );
    void FC_FUNC_(iaga_get_eigenvalues, IAGA_GET_EIGENVALUES)( double*, int* );
  }
  //! \endcond
#endif


namespace Pescan
{
# ifdef _DIRECTIAGA
    void Interface :: set_mpi( boost::mpi::communicator* _c )
    {
      MPI_COMMDEC::set_mpi( _c );
      MPI_Comm __commC = (MPI_Comm) ( MPI_COMM ) ;
      MPI_Fint __commF = MPI_Comm_c2f( __commC );
      FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );
    }
# endif

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
    const t_Path newpath(  opt::InitialPath::path() / dirname );
    __ASSERT( newpath.empty(), "Path is empty.\n" )
    __DIAGA( for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i ) )
    {
      __DIAGA
      ( 
        MPI_COMM.barrier(); 
        if( MPI_COMM.rank() != i ) continue;
      )
      if( not bfs::exists( newpath ) ) bfs::create_directory( newpath );
    }
  }
  void Interface :: destroy_directory()
  {
    namespace bfs = boost::filesystem;
    const t_Path newpath(  opt::InitialPath::path() / dirname );
    __DIAGA( for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i ) )
    {
      __DIAGA
      ( 
        MPI_COMM.barrier(); 
        if( MPI_COMM.rank() == i ) continue;
      )
      if( bfs::exists( newpath ) ) bfs::remove_all( newpath );
    }
  }
  void Interface :: create_potential()
  {
    namespace bfs = boost::filesystem;
    write_genpot_input();
    
    const t_Path newatom( opt::InitialPath::path()/dirname/atom_input.filename() );
    __ASSERT( not bfs::exists( atom_input ), atom_input << " does not exist.\n" );
    if( newatom != atom_input )
    { 
      if( bfs::exists( newatom ) ) bfs::remove( newatom );
      bfs::copy_file( atom_input, newatom );
    }
    __ASSERT( not bfs::exists( newatom ), newatom << " does not exist.\n" );

#   ifndef _NOLAUNCH
      __IIAGA
      ( 
        const t_Path newgenpot( opt::InitialPath::path()/dirname/genpot.launch.filename() );
        if( not bfs::exists( newgenpot ) ) bfs::copy_file( genpot.launch, newgenpot ); 
      )
#   endif
    __DIAGA( for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i ) )
    {
      __DIAGA
      (
        MPI_COMM.barrier();
        if( MPI_COMM.rank() != i ) continue;
      )
      foreach( const t_Path& ppath,  genpot.pseudos )
      {
        const t_Path newstr( opt::InitialPath::path()/dirname/ppath.filename() );
        if( not bfs::exists( newstr ) ) bfs::copy_file( ppath, newstr );
      }
    }

    
    chdir( ( opt::InitialPath::path() / dirname ).string().c_str() );
#   ifndef _NOLAUNCH
      __IIAGA(  system( "./" + genpot.launch.filename().string() ) );
      __DIAGA( FC_FUNC_(getvlarg, GETVLARG)(); )
#   endif
    chdir( opt::InitialPath::path().string().c_str() );
  }
  types::t_real Interface :: launch_pescan()
  {
    namespace bfs = boost::filesystem;
    write_escan_input();

    __DIAGA( for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i ) )
    {
      __DIAGA
      ( 
        MPI_COMM.barrier(); 
        if( MPI_COMM.rank() != i ) continue;
      )
      const t_Path newmaskr( opt::InitialPath::path()/dirname/maskr.filename() );
      if( not bfs::exists( newmaskr ) ) bfs::copy_file( maskr, newmaskr );
    }
#   ifndef _NOLAUNCH
      __IIAGA
      (
        const t_Path newescan( opt::InitialPath::path()/dirname/escan.launch.filename() );
        if( not bfs::exists( newescan ) ) bfs::copy_file( escan.launch, newescan );
      )
#   endif
    foreach( const SpinOrbit &sp, escan.spinorbit )
      __DIAGA( for( types::t_int i( 0 ); i < MPI_COMM.size(); ++i ) )
      {
        __DIAGA
        ( 
          MPI_COMM.barrier(); 
          if( MPI_COMM.rank() != i ) continue;
        )
        const t_Path newso( opt::InitialPath::path()/dirname/ sp.filename.filename() );
        if( not bfs::exists( newso ) ) bfs::copy_file( sp.filename, newso );
      }
         
    chdir( ( opt::InitialPath::path() / dirname ).string().c_str() );
#   ifndef _NOLAUNCH
      __IIAGA(  system( "./" + escan.launch.filename().string() ) );
      __DIAGA( FC_FUNC_(iaga_call_escan, IAGA_CALL_ESCAN)( &escan.nbstates ); )
#   endif
    chdir( opt::InitialPath::path().string().c_str() );
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

    __DOASSERT( not genpot.Load( _node ), "Could not load genpot input.\n" );

    do_destroy_dir = true;
    if( _node.Attribute("keepdirs") ) do_destroy_dir = false;

    child = _node.FirstChildElement("Maskr");
    if( child and child->Attribute("filename") )
      maskr = child->Attribute("filename");

    __DOASSERT( not escan.Load( _node ), "Could not load escan input.\n" );

    __TRYCODE( check_existence();, "Some files and/or programs are missing.\n" )

    return true;
  }

  void Interface::write_genpot_input()
  {
    std::ofstream file;
    const t_Path orig(   opt::InitialPath::path() 
                       / dirname / __DIAGASUFFIX( genpot.filename ) );
    file.open( orig.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 

    __DOASSERT( file.bad() or ( not file.is_open() ),
                "Could not open file " << orig << " for writing.\nAborting.\n" )

    __DOASSERT( not genpot.check(), "Input to genpot is not coherent.\n" )

    file << Interface :: t_Path(atom_input.stem()).filename() << "\n"
         << genpot;

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
        __DIAGA( __DIAGASUFFIX( t_Path("escan_input") ) )
    );
    file.open( orig.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 

    __DOASSERT( file.bad() or ( not file.is_open() ),
                "Could not open file " << orig << " for writing.\nAborting.\n" )

    file << "1 " << __DIAGASUFFIX( genpot.output ) << "\n" 
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
      file << "\n" << "9 " << escan.wavefunction_in << "\n";
    }
    else  file << "7 0\n8 0\n9 dummy\n";
    file << "10 " << escan.rspace_output << " "
                  << " 1 1 1 " << escan.rspace_wfn << "\n";

    if ( atat::norm2( escan.kpoint ) < types::tolerance ) file << "11 0 0 0 0 0\n";
    else
    {
      atat::rVector3d k = escan.kpoint;
      file << "11 1 " << k << " " 
                      << std::setw(12) << std::setprecision(8) 
                      << escan.scale / Physics::a0("A") <<  "\n";
    }
    file << "12 " << escan.potential << "\n"
         << "13 " << atom_input.filename() << "\n"
         << "14 " << escan.rcut << "\n"
         << "15 " << escan.spinorbit.size() << "\n";
    std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
    std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
    for(types::t_unsigned i=16; i_so != i_so_end; ++i_so, ++i )
      file << i << " " << i_so->filename.filename() << " " 
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
        eigenvalues.resize( escan.nbstates );
        std::fill( eigenvalues.begin(), eigenvalues.end(), 0 );
        FC_FUNC_
        ( 
          iaga_get_eigenvalues, IAGA_GET_EIGENVALUES
        )(&eigenvalues[0],&escan.nbstates);
        boost::mpi::broadcast( MPI_COMM, eigenvalues, 0 );
      )
     
      __IIAGA( __DOASSERT( u != escan.nbstates,
                              "Found " << u << " eigenvalues in " << name
                           << " where " << escan.nbstates
                           << " were expected.\n" ) )
#   endif
    return true;
  }
               
  bool Interface :: GenPot :: check()
  { 
    if( multiple_cell == t_MeshTuple(0,0,0) ) 
    {
      small_box = multiple_cell;
      multiple_cell = mesh;
    }
    else if(    boost::tuples::get<0>(mesh)
                  % boost::tuples::get<0>(multiple_cell) != 0 
             or boost::tuples::get<1>(mesh)
                  % boost::tuples::get<1>(multiple_cell) != 0  
             or boost::tuples::get<2>(mesh)
                  % boost::tuples::get<2>(multiple_cell) != 0  ) return false;
    return pseudos.size() >= 2; 
  }

  bool Interface :: GenPot :: Load( const TiXmlElement &_node )
  {
    try 
    {
      const TiXmlElement *child = _node.FirstChildElement("GenPot");
      __DOASSERT( not child, "No <GenPot> tag found on input.\nAborting\n" )
      __DOASSERT( not child->Attribute("x"), "x attributes are missing in <GenPot>.\n")
      __DOASSERT( not child->Attribute("y"), "y attributes are missing in <GenPot>.\n")
      __DOASSERT( not child->Attribute("z"), "z attributes are missing in <GenPot>.\n")
      __DOASSERT( not child->Attribute("cutoff"),
                  "cutoff attributes are missing in <GenPot>.\n")
      __DOASSERT( not child->FirstChildElement("Pseudo"),
                  "Could not find <Pseudo/> in <GenPot>\nAborting\n")
     
      boost::tuples::get<0>(mesh)
        = boost::lexical_cast<types::t_real>( child->Attribute("x") );
      boost::tuples::get<1>(mesh)
          = boost::lexical_cast<types::t_real>( child->Attribute("y") );
      boost::tuples::get<2>(mesh)
          = boost::lexical_cast<types::t_real>( child->Attribute("z") );
      multiple_cell = mesh;
      if( not ( child->Attribute("mx") or child->Attribute("my") or child->Attribute("mz") ) )
      {
        boost::tuples::get<0>(multiple_cell)
          = boost::lexical_cast<types::t_real>( child->Attribute("mx") );
        boost::tuples::get<1>(multiple_cell)
          = boost::lexical_cast<types::t_real>( child->Attribute("my") );
        boost::tuples::get<2>(multiple_cell)
          = boost::lexical_cast<types::t_real>( child->Attribute("mz") );
      }
      small_box = t_MeshTuple(0,0,0);
      if( not ( child->Attribute("shmx") or child->Attribute("shmy") or child->Attribute("shmz") ) )
      {
        boost::tuples::get<0>(small_box)
          = boost::lexical_cast<types::t_real>( child->Attribute("shmx") );
        boost::tuples::get<1>(small_box)
          = boost::lexical_cast<types::t_real>( child->Attribute("shmy") );
        boost::tuples::get<2>(small_box)
          = boost::lexical_cast<types::t_real>( child->Attribute("shmz") );
      }
      child->Attribute("cutoff", &cutoff);
      __IIAGA
      (
        if ( child->Attribute("launch") )
          launch = child->Attribute("launch");
      )
      child = child->FirstChildElement("Pseudo");
      for(; child; child = child->NextSiblingElement() )
        if ( child->Attribute("filename") )
          pseudos.push_back( child->Attribute("filename") );
      __DOASSERT( not check(),
                  "Insufficient or Incorrect Genpot input.\nAborting\n" )
    }
    catch( std::exception &_e )
    {
      std::cerr << "Error while loading Genpot from XML.\n" << _e.what() << "\n";
      return false;
    }
    return true;
  } 

  std::ostream& operator<<( std::ostream& _stream, const Interface::GenPot &_g )
  {
    _stream << _g.mesh << "\n" 
            << _g.multiple_cell << "\n" 
            << _g.small_box << "\n"
            << _g.cutoff << "\n"
            << _g.pseudos.size() << std::endl;
    foreach( const Interface :: t_Path& ppath, _g.pseudos )
      _stream << ppath.filename() << "\n";
    return _stream;
  }

  bool Interface :: Escan :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *child = _node.FirstChildElement("Hamiltonian");
    try
    {
      __DOASSERT( not child, "Could not find <Hamiltonian> tag on input.\n")
      __DOASSERT( not child->Attribute("kinscal"),
                  "kinscal attribute is missing in <Hamiltonian> tag.\n")
      __DOASSERT( not child->Attribute("realcutoff"),
                  "realcutoff attribute is missing in <Hamiltonian> tag.\n")
      __DOASSERT( not child->Attribute("potential"),
                  "potential attribute is missing in <Hamiltonian> tag.\n")
    }
    catch( std::exception& _e )
    {
      std::cerr << "Escan input is incomplete.\n" << _e.what();
      return false;
    }
    __IIAGA
    (
      if ( child->Attribute("launch") )
        launch = child->Attribute("launch");
    )
    child = _node.FirstChildElement("Wavefunctions");
    if( child and child->Attribute("in") )
      wavefunction_in = child->Attribute("in");
    if( child and child->Attribute("out") )
      wavefunction_out = child->Attribute("out");
    if( _node.Attribute("method") )
    {
      const std::string method_string( child->Attribute("method") );
      method = ALL_ELECTRON;
      if(    method_string.find("folded") != std::string::npos
          or method_string.find("1") != std::string::npos )
        method = FOLDED_SPECTRUM;
    }
    child = _node.FirstChildElement("Reference");
    if( child and child->Attribute("value") )
      child->Attribute("value", &Eref);



    child->Attribute("kinscal", &kinscal);
    if( child->Attribute("smooth") )
      child->Attribute("smooth", &smooth);
    if( child->Attribute("potential") )
    {
      const std::string str( child->Attribute("potential") );
      if(    str.find( "local" ) != std::string::npos 
          or str.find( "1" ) != std::string::npos  )
        potential = Escan::LOCAL; 
      else if(    str.find( "nonlocal" ) != std::string::npos 
               or str.find( "non-local" ) != std::string::npos  
               or str.find( "2" ) != std::string::npos  )
        potential = Escan::NONLOCAL; 
      else if(    str.find( "spinorbit" ) != std::string::npos 
               or str.find( "spin-orbit" ) != std::string::npos  
               or str.find( "3" ) != std::string::npos  )
        potential = Escan::SPINORBIT; 
    }
    if( child->Attribute("realcutoff") )
      child->Attribute("realcutoff", &rcut);
    if( child->Attribute("nbstates") )
      child->Attribute("nbstates", &nbstates);
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
        spinorbit.push_back(so);
      }
    child = _node.FirstChildElement("Minimizer");
    if( not child ) return true;

    if( child->Attribute("niter") )
      child->Attribute("niter", &niter);
    if( child->Attribute("nlines") )
      child->Attribute("nlines", &nlines);
    if( child->Attribute("tolerance") )
      child->Attribute("tolerance", &tolerance);

    return true;
  }
}
