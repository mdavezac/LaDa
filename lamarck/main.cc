//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>
#include <sstream>
#include <string>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/program_options.hpp>

#include <print/manip.h>
#include <mpi/mpi_object.h>

#include "functional_builder.h"

#include <revision.h>
#define __PROGNAME__ "lamarck"

#ifdef _LADADEBUG
#include <print/stdout.h>
#define OUTPUT Print::out
#define ENDLINE Print::endl
#else
#define OUTPUT std::cout
#define ENDLINE "\n"
#endif
#ifndef TIXML_USE_STL
#error not using TIXML_USE_STL
#endif

#ifdef _MPI
#include <boost/mpi/environment.hpp>
#endif

#include "constituent_strain.h"
#include "harmonic.h"

#if defined(_CUBIC_CE_)
typedef Ising_CE::ConstituentStrain::Harmonic::Cubic t_Harmonic;
#elif defined( _TETRAGONAL_CE_ )
typedef Ising_CE::ConstituentStrain::Harmonic::Tetragonal t_Harmonic;
#else
#error Please specify _CUBIC_CE_ or _TETRAGONAL_CE_
#endif
typedef VA_CE::Builder<t_Harmonic> t_Builder;

int main(int argc, char *argv[]) 
{
  try
  {
    __MPICODE(
      boost::mpi::environment env(argc, argv); 
      boost::mpi::communicator world;
      ::mpi::main = &world;
      Print::out.init( "out" );
      Print::out.doprint( true );
    )

    std::string filename("input.xml");
    namespace po = boost::program_options;

    po::options_description generic("Generic Options");
    generic.add_options()
           ("help,h", "produces this help message.")
           ("version,v", "prints version string.");
    po::options_description specific("GA Options");
    specific.add_options()
        ("input,i", po::value<std::string>()->default_value("input.xml"), "input filename." )
        ("reruns,r", po::value<unsigned>()->default_value(1),
                     "number of times to run the algorithm.\n"
                     "Is equivalent to manually re-launching the program.\n");
 
    po::options_description all;
    all.add(generic).add(specific);
 
    po::positional_options_description p;
    p.add("input", 1);
 
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(p).run(), vm);
    po::notify(vm);
 
    if ( vm.count("help") )
    {
      __ROOTCODE( (*::mpi::main),
        OUTPUT << "Usage: " << argv[0] << " [options] file.xml\n"
                  << "  file.xml is an optional filename for XML input.\n"
                  << "  Default input is input.xml.\n\n"
                  << all << "\n"; 
      )
      return 1;
    }
    if ( vm.count("version") )
    {
      __ROOTCODE( (*::mpi::main),
        OUTPUT << "\n" << __PROGNAME__
                  << " from the " << PACKAGE_STRING << " package\n"
                  << "Subversion Revision: " << SVN::Revision << "\n\n"; 
      )
      return 1;
    }
    if ( vm.count("input") )
    {
      filename = vm["input"].as< std::string >();
      if( filename.compare("input.xml") )
        __ROOTCODE( (*::mpi::main), OUTPUT << "Input: " << filename << ".\n";)
    }

    TiXmlElement *child;
    atat::rVector3d vec;
    Ising_CE::Lattice lattice;
 
    
    __MPICODE(
      Print :: out << "Before broadcasting" << Print::endl;
      boost::mpi::broadcast( *::mpi::main, filename, 0 );
      Print :: out << "After broadcasting" << Print::endl;
      std::string input;
    )
 
    __ROOTCODE( (*::mpi::main),
      TiXmlDocument doc( filename );
      __DOASSERT( not doc.LoadFile(), 
                    "error while opening input file "
                 << filename << "\n" << doc.ErrorDesc()  )
      __MPICODE( 
        std::ostringstream stream;
        stream << *doc.RootElement();
        input = stream.str();
      )
    )
    
    __MPICODE(
        boost::mpi::broadcast( *::mpi::main, input, 0 ); 
        TiXmlDocument doc;
        doc.Parse( input.c_str() );
    )
    TiXmlHandle handle( &doc );
 
    // loads lattice
    child = handle.FirstChild( "Job" ).FirstChild( "Lattice" ).Element();
    __DOASSERT( not child, "Could not find Lattice in input." )
    __DOASSERT( not lattice.Load(*child), "Error while reading Lattice from input.")
 
    // loads lamarck functional
    child = handle.FirstChild( "Job" ).Element();
    __DOASSERT( not child, "Could not find Functional in input." )

    t_Builder ce;
    __DOASSERT( not ce.Load(*child), "Error while reading Functional from input." )
    ce.add_equivalent_clusters();
 
    __MPIROOT( (*::mpi::main), OUTPUT << "Nb procs: " << ::mpi::main->size() << ENDLINE; )
    // do structure
    child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
    for (; child; child = child->NextSiblingElement("Structure") )
    {
      Ising_CE::Structure structure;
      Ising_CE :: Structure :: lattice = &lattice;
      __DOASSERT( not structure.Load(*child), "Error while reading Structure from input." )
 
      t_Builder::t_VA_Functional functional;
      ce.generate_functional(structure, &functional);
      __MPICODE( functional.get_functional1()->set_mpi( ::mpi::main ); )
      __MPICODE( functional.get_functional2()->set_mpi( ::mpi::main ); )
    
      functional.resize( structure.atoms.size() );
      std::transform( structure.atoms.begin(), structure.atoms.end(), functional.begin(),
                      boost::lambda::bind( &Ising_CE::Structure::t_Atom::type,
                                           boost::lambda::_1 ) );
    
      OUTPUT << "Energy: " << functional.evaluate() // << "\n";
             << "  Concentration: " << structure.get_concentration() << "\n" ; //<< ENDLINE;
//   Ising_CE::Fourier( structure.atoms.begin(), structure.atoms.end(),
//                      structure.k_vecs.begin(), structure.k_vecs.end() );
//   structure.print_out( std::cout );

//   delete functional.get_functional1();
//   delete functional.get_functional2();
    }
  }
  catch ( boost::program_options::invalid_command_line_syntax &_b)
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Something wrong with the command-line input.\n"
              << _b.what() << "\n";
    return 0;
  }
  catch ( boost::program_options::invalid_option_value &_i )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Argument of option in command-line is invalid.\n"
              << _i.what() << "\n";
    return 0;
  }
  catch ( boost::program_options::unknown_option &_u)
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Unknown option in command-line.\n"
              << _u.what() << "\n";
    return 0;
  }
  catch (  boost::program_options::too_many_positional_options_error &_e )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ << "\n"
              << "Too many arguments in command-line.\n"
              << _e.what() << "\n";
    return 0;
  }
  catch ( std::exception &e )
  {
    std::ostringstream sstr;
    __MPICODE( sstr << "\nProcessor " << mpi::main->rank() + 1
                    << " of " << mpi::main->size()
                    << " says:\n"; )

    sstr << "Caught error while running " << __PROGNAME__ 
         << "\n" << e.what() << "\n";

    std::string message = sstr.str();
    __MPICODE( 
      boost::mpi::reduce( *::mpi::main, sstr.str(), message, std::plus<std::string>(), 0 );
    )
    __NOTMPIROOT( (*::mpi::main), return 0; )
    std::cerr << message << std::endl;
    __DODEBUGCODE( Print::out << message << Print::endl; )
    return 0;
  }
  return 1;
}

