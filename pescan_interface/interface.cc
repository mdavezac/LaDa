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

#include <boost/xpressive/regex_primitives.hpp>
#include <boost/xpressive/regex_algorithms.hpp>
#include <boost/xpressive/regex_compiler.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <opt/debug.h>
#include <opt/initial_path.h>
#include <opt/tuple_io.h>
#include <opt/tinyxml.h>
#include <crystal/lattice.h>
#include <print/manip.h>
#include <mpi/macros.h>

#include "interface.h"


#ifdef _DIRECTIAGA
  // declares fortran interface
  //! \cond
  extern "C"
  {
    void FC_FUNC_(getvlarg, GETVLARG)();
    void FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( MPI_Fint * );
    void FC_FUNC_(iaga_call_escan, IAGA_CALL_ESCAN)( int*, int* );
    void FC_FUNC_(iaga_get_eigenvalues, IAGA_GET_EIGENVALUES)( double*, int* );
  }
  //! \endcond
#endif


namespace LaDa
{
  namespace Pescan
  {
#   ifdef _DIRECTIAGA
      void Interface :: set_mpi( boost::mpi::communicator* _c )
      {
        MPI_COMMDEC::set_mpi( _c );
        MPI_Comm __commC = (MPI_Comm) ( MPI_COMM ) ;
        MPI_Fint __commF = MPI_Comm_c2f( __commC );
        FC_FUNC_(iaga_set_mpi, IAGA_SET_MPI)( &__commF );
      }
#   endif

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
        int verb( verbose ? 1: 0);
        __DIAGA( FC_FUNC_(iaga_call_escan, IAGA_CALL_ESCAN)( &escan.nbstates, &verb ); )
  #   endif
      chdir( opt::InitialPath::path().string().c_str() );
      return 0.0;
    }

 
    bool Interface::Load( const TiXmlElement &_node )
    {
      namespace bfs = boost::filesystem;
      const TiXmlElement *parent = opt::find_node( _node, "Functional", "type", "escan" );

      if( not parent )
      {
        std::cerr << "Could not find an <Functional type=\"escan\"> tag in input file" 
                  << std::endl;
        return false;
      }
      if( not parent->Attribute( "filename" ) ) return Load_( *parent );

      const bfs::path path( Print::reformat_home( parent->Attribute( "filename" ) ) );
      __DOASSERT( not bfs::exists( path ), path.string() + " does not exist.\n" )
      TiXmlDocument doc;
      opt::read_xmlfile( path, doc );
      __DOASSERT( not doc.FirstChild( "Job" ),
                  "Root tag <Job> does not exist in " + path.string() + ".\n" )
      parent = opt::find_node( *doc.FirstChildElement( "Job" ), "Functional", "type", "escan" );
      
      if( parent ) return Load_( *parent );

      std::cerr << "Could not find an <Functional type=\"escan\"> tag in input file" 
                << std::endl;
      return false;
    }

    bool Interface :: Load_ (const TiXmlElement &_node )
    {
      const TiXmlElement *child;

      // lester needs to load "module fftw" when in interactive multi-proc mode
  #   ifndef _DIRECTIAGA
        if ( _node.Attribute("system") )
        {
          std::string str = _node.Attribute("system");
          system( str.c_str() );
        }
  #   endif
      if( _node.Attribute("method") )
      {
        const std::string method_string( _node.Attribute("method") );
        escan.method = ALL_ELECTRON;
        if(    method_string.find("folded") != std::string::npos
            or method_string.find("1") != std::string::npos )
          escan.method = FOLDED_SPECTRUM;
      }
      verbose = true;
      if( _node.Attribute("verbose") )
      {
        namespace bx = boost::xpressive;
        std::string const verb_string = _node.Attribute("verbose");
        bx::sregex false_ = bx::icase( bx::as_xpr('f') >> !bx::as_xpr("alse") ) | '1' ;
        bx::sregex true_ = bx::icase( bx::as_xpr('t') >> !bx::as_xpr("rue") ) | '0' ;
        if( bx::regex_match( verb_string, false_ ) ) verbose = false;
        else if( bx::regex_match( verb_string, true_ ) ) verbose = true;
        else __THROW_ERROR( "Uknown argument to verbose in Pescan: " + verb_string + "\n" )
      }

      __DOASSERT( not genpot.Load( _node ), "Could not load genpot input.\n" );

      do_destroy_dir = true;
      if( _node.Attribute("keepdirs") ) do_destroy_dir = false;

      child = _node.FirstChildElement("Maskr");
      if( child and child->Attribute("filename") )
        maskr = Print::reformat_home( child->Attribute("filename") );

      __DOASSERT( not escan.load_hamiltonian( _node ), "Could not load escan input.\n" );

      if( escan.load_reference( _node ) ) escan.method = FOLDED_SPECTRUM;
      escan.load_wavefunctions( _node );
      escan.load_minimizer( _node );

      __TRYCODE( check_existence();, "Some files and/or programs are missing.\n" )

      return true;
    }

    void Interface::write_genpot_input()
    {
      namespace bt = boost::tuples;
      std::ofstream file;
      const t_Path orig(   opt::InitialPath::path() 
                         / dirname / __DIAGASUFFIX( genpot.filename ) );
      file.open( orig.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 

      __DOASSERT( file.bad() or ( not file.is_open() ),
                  "Could not open file " << orig << " for writing.\nAborting.\n" )

      __TRYCODE( genpot.check();, "Input to genpot is not coherent.\n" )
      __DIAGA
      (
        __DOASSERT
        (
             bt::get<0>( genpot.mesh ) * bt::get<1>(  genpot.mesh ) * bt::get<2>( genpot.mesh )
           % MPI_COMM.size() != 0,
              "MPI pool size is inadequate: pool \% n1 * n2 * n3 != 0\n"
           << bt::get<0>( genpot.mesh ) << "*" << bt::get<1>( genpot.mesh ) 
           << "*" << bt::get<2>( genpot.mesh ) << " \% " << MPI_COMM.size()
           << " = " << (   bt::get<0>( genpot.mesh ) * bt::get<1>( genpot.mesh )
                         * bt::get<2>( genpot.mesh ) % MPI_COMM.size() != 0 )
        )
      )

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
        / __IIAGA( "escan.input" )
          __DIAGA( __DIAGASUFFIX( t_Path("escan_input") ) )
      );
      file.open( orig.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 

      __DOASSERT( file.bad() or ( not file.is_open() ),
                  "Could not open file " << orig << " for writing.\nAborting.\n" )

      file << "1 " << __DIAGASUFFIX( genpot.output ) << "  # generated atomic potential\n" 
           << "2 " << escan.wavefunction_out << " # output G-space wavefunctions\n"
           << "3 " << escan.method << "  # escan method: "
           << ( escan.method == FOLDED_SPECTRUM  ?
                   " 1 = folded-spectrum\n": " 2 = diagonalization\n" )
           << "4 " << escan.Eref << " " << genpot.cutoff << " "
                   << escan.smooth << " " << escan.kinscal 
                   << " # Eref, cutoff, smooth, kinetic energy scaling\n"
           << "5 " << escan.nbstates << " # number of states\n"
           << "6 " << escan.niter << " " << escan.nlines
                   << " " << escan.tolerance << " # nb iterations, nlines, tolerance \n";
      if( escan.read_in.size() )
      {
        file << "7 " << escan.read_in.size() << " # nb of wfns to read from input\n"
             << "8 ";
        std::vector<types::t_unsigned> :: const_iterator i_r = escan.read_in.begin();
        std::vector<types::t_unsigned> :: const_iterator i_r_end = escan.read_in.end();
        file << *i_r;
        for(--i_r_end; i_r != i_r_end; ++i_r)
          file << ", " << *i_r;
        file << " # wfn indices \n" << "9 " << escan.wavefunction_in << " # wfn input filename\n";
      }
      else  file << "7 0 # nb of wfns read on input \n8 0 # unused \n9 dummy # unused \n";
      file << "10 " << escan.rspace_output << " "
                    << " 1 1 1 " << escan.rspace_wfn << " # real-space output \n";

      if ( math::is_zero( escan.kpoint.squareNorm() ) ) file << "11 0 0 0 0 0 # kpoint=Gamma\n";
      else
      {
        math::rVector3d k = escan.kpoint;
        file << "11 1 " << k << " " 
                        << std::setw(12) << std::setprecision(8) 
                        << escan.scale / Physics::a0("A") <<  " # kpoint\n";
      }
      file << "12 " << escan.potential << " # potential:";
      switch( escan.potential )
      {
        case Escan::NOPOT: __THROW_ERROR( "unknown potential." ); break;
        case Escan::LOCAL: file << " local potential\n"; break;
        case Escan::NONLOCAL: file << " non-local potential\n"; break;
        case Escan::SPINORBIT: file << " non-local potential with spin-orbit coupling\n"; break;
      }
      file << "13 " << atom_input.filename() << " # atomic configuration filename\n"
           << "14 " << escan.rcut << " # real-space projector cut-off\n"
           << "15 " << escan.spinorbit.size() << " # number of spin-orbit potentials\n";
      std::vector<SpinOrbit> :: const_iterator i_so = escan.spinorbit.begin();
      std::vector<SpinOrbit> :: const_iterator i_so_end = escan.spinorbit.end();
      for(types::t_unsigned i=16; i_so != i_so_end; ++i_so, ++i )
        file << i << " " << i_so->filename.filename() << " " 
             << i_so->izz << " " 
             << i_so->s << " " << i_so->p << " " << i_so->d << " " 
             << i_so->pnl << " " << i_so->dnl 
             << " # filename, type, s, p, d, p-nonlocal, d-nonlocal.\n";
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
                 
    void Interface :: GenPot :: check()
    { 
      namespace bt = boost::tuples;
      __DOASSERT
      (     
           bt::get<0>( mesh ) == 0        
        or bt::get<1>( mesh ) == 0      
        or bt::get<2>( mesh ) == 0,
        "Requested mesh has a null component.\n"
      )
      __DOASSERT
      (     
           bt::get<0>( multiple_cell ) == 0         
        or bt::get<1>( multiple_cell ) == 0      
        or bt::get<2>( multiple_cell ) == 0,
        "Requested multiple_cell has a null component.\n"
      )
      __DOASSERT
      (     
           bt::get<0>( mesh ) % bt::get<0>( multiple_cell ) != 0         
        or bt::get<1>( mesh ) % bt::get<1>( multiple_cell ) != 0      
        or bt::get<2>( mesh ) % bt::get<2>( multiple_cell ) != 0,
           "Requested multiple cell does not evenly divide requested mesh.\n"
        << mesh << " % " << multiple_cell << "\n"
      )
      __DOASSERT
      (     
           bt::get<0>( mesh ) < bt::get<0>( multiple_cell )         
        or bt::get<1>( mesh ) < bt::get<1>( multiple_cell )      
        or bt::get<2>( mesh ) < bt::get<2>( multiple_cell ),
        "Requested multiple cell parameters are larger than the requested mesh.\n"
      )
      __DOASSERT
      (     
           bt::get<0>( multiple_cell ) < bt::get<0>( small_box )         
        or bt::get<1>( multiple_cell ) < bt::get<1>( small_box )      
        or bt::get<2>( multiple_cell ) < bt::get<2>( small_box ),
        "Requested multiple cell parameters are smaller than the small box parameters.\n"
      )
      __DOASSERT( pseudos.size() < Crystal::nb_species( *Crystal::Structure::lattice ),
                  "Fewer pseudos than species in lattice.\n" )
    }

    bool Interface :: GenPot :: Load( const TiXmlElement &_node )
    {
      namespace bt = boost::tuples;
      try 
      {
        const TiXmlElement *child = _node.FirstChildElement("GenPot");
        __DOASSERT( not child, "No <GenPot> tag found on input.\nAborting\n" )
        if( not child->Attribute("mesh") )
        {
          __DOASSERT( not child->Attribute("x"), "x attributes are missing in <GenPot>.\n")
          __DOASSERT( not child->Attribute("y"), "y attributes are missing in <GenPot>.\n")
          __DOASSERT( not child->Attribute("z"), "z attributes are missing in <GenPot>.\n")
        }
        __DOASSERT( not child->Attribute("cutoff"),
                    "cutoff attributes are missing in <GenPot>.\n")
        __DOASSERT( not child->FirstChildElement("Pseudo"),
                    "Could not find <Pseudo/> in <GenPot>\nAborting\n")
       
        if( child->Attribute("mesh") )
        {
          __DOASSERT( not opt::tuples::read( child->Attribute("mesh"), mesh ),
                      "Could not parse mesh attribute.\n" )
        }
        else 
        {
          bt::get<0>(mesh) = boost::lexical_cast<types::t_int>( child->Attribute("x") );
          bt::get<1>(mesh) = boost::lexical_cast<types::t_int>( child->Attribute("y") );
          bt::get<2>(mesh) = boost::lexical_cast<types::t_int>( child->Attribute("z") );
        }
        multiple_cell = mesh;
        if( child->Attribute("multcell") )
        {
          __DOASSERT( not opt::tuples::read( child->Attribute("multcell"), multiple_cell ),
                      "Could not parse multcell attribute.\n" )
        }
        small_box = multiple_cell == mesh ? t_MeshTuple(0,0,0): multiple_cell;
        if( child->Attribute("smallbox") )
        {
          __DOASSERT( not opt::tuples::read( child->Attribute("smallbox"), small_box ),
                      "Could not parse smallbox attribute.\n" )
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
            pseudos.push_back( Print::reformat_home( child->Attribute("filename") ) );
        check();
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
      _stream << boost::tuples::set_open(' ')
              << boost::tuples::set_close(' ')
              << boost::tuples::set_delimiter(' ')
              << _g.mesh << "\n" 
              << boost::tuples::set_open(' ')
              << boost::tuples::set_close(' ')
              << boost::tuples::set_delimiter(' ')
              << _g.multiple_cell << "\n" 
              << boost::tuples::set_open(' ')
              << boost::tuples::set_close(' ')
              << boost::tuples::set_delimiter(' ')
              << _g.small_box << "\n"
              << _g.cutoff << "\n"
              << _g.pseudos.size() << std::endl;
      foreach( const Interface :: t_Path& ppath, _g.pseudos )
        _stream << ppath.filename() << "\n";
      return _stream;
    }

    bool Interface :: Escan :: load_hamiltonian( const TiXmlElement &_node )
    {
      // Reads attribute of escan functional tag.
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
      if( child->Attribute("kpoint") )
      {
        boost::tuples::tuple<types::t_real, types::t_real, types::t_real> kp;
        if( not opt::tuples::read( child->Attribute("kpoint"), kp ) )
          std::cerr << "Could not parse kpoint attribute in Hamiltonian tag.\n";
        else
        {
          kpoint[0] = boost::tuples::get<0>(kp);
          kpoint[1] = boost::tuples::get<1>(kp);
          kpoint[2] = boost::tuples::get<2>(kp);
        }
      }
      child = child->FirstChildElement("SpinOrbit");
      for(; child; child = child->NextSiblingElement("SpinOrbit") )
        if ( child->Attribute("filename") and child->Attribute("izz") )
        {
          SpinOrbit so;
          so.filename = Print::reformat_home( child->Attribute("filename") );
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

      return true;
    }

    bool Interface :: Escan :: load_wavefunctions( const TiXmlElement &_node )
    {
      // Reads attribute of escan functional tag.

      const TiXmlElement *child = _node.FirstChildElement("Wavefunctions");
      if( not child ) return false;
      if( child->Attribute("in") )
        wavefunction_in = child->Attribute("in");
      if( child->Attribute("out") )
        wavefunction_out = child->Attribute("out");

      return true;
    }

    bool Interface :: Escan :: load_minimizer( const TiXmlElement &_node )
    {
      // Reads attribute of escan functional tag.
      const TiXmlElement *child = _node.FirstChildElement("Minimizer");
      if( not child ) return true;

      if( child->Attribute("niter") )
        niter = boost::lexical_cast<types::t_int>( child->Attribute("niter") );
      if( child->Attribute("nlines") )
        nlines = boost::lexical_cast<types::t_int>( child->Attribute("nlines") );
      if( child->Attribute("tolerance") )
        tolerance = boost::lexical_cast<types::t_real>( child->Attribute("tolerance") );

      return true;
    }

    bool Interface :: Escan :: load_reference( const TiXmlElement &_node )
    {
      // Reads attribute of escan functional tag.
      const TiXmlElement *child = _node.FirstChildElement("Reference");
      if( not child ) return false;
      if( child->Attribute("value") )
        Eref = boost::lexical_cast< types::t_real >( child->Attribute("value") );
      return true;
    }
    
    size_t nb_valence_states( const Crystal::Structure &_str ) 
    {
      Crystal::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
      Crystal::Structure::t_Atoms::const_iterator i_atom_end = _str.atoms.end();
      types::t_unsigned bgstates = 0;
      for(; i_atom != i_atom_end; ++i_atom)
      {
        Crystal::StrAtom atom; 
        _str.lattice->convert_Atom_to_StrAtom( *i_atom, atom );
        bgstates += Physics::Atomic::Charge( atom.type );
      }
      return bgstates;
    }
  }
} // namespace LaDa
