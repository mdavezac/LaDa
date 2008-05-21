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

int main(int argc, char *argv[]) 
{
  try
  {
    mpi::environment env(argc, argv);

    namespace po = boost::program_options;
    std::string filename("input.xml");

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
      __ROOTCODE(
        std::cout << "Usage: " << argv[0] << " [options] file.xml\n"
                  << "  file.xml is an optional filename for XML input.\n"
                  << "  Default input is input.xml.\n\n"
                  << all << "\n"; 
      )
      return 1;
    }
    if ( vm.count("version") )
    {
      __ROOTCODE( 
        std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package\n"
                  << "Subversion Revision: " << SVN::Revision << "\n\n"; 
      )
      return 1;
    }
    if ( vm.count("input") )
    {
      filename = vm["input"].as< std::string >();
      __ROOTCODE(std::cout << "Input: " << filename << ".\n";)
    }

    TiXmlElement *child;
    atat::rVector3d vec;
    Ising_CE::Lattice lattice;
 
    
    boost::mpi::broadcast( *::mpi::main, filename, 0 );
    TiXmlDocument doc( filename.c_str() );
    
    if  ( !doc.LoadFile() )
    {
      std::cerr << "error while opening input file " << filename << std::endl
                << doc.ErrorDesc() << std::endl; 
      return false;
    }
 
    TiXmlHandle handle( &doc );
 
    // loads lattice
    child = handle.FirstChild( "Job" ).FirstChild( "Lattice" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find Lattice in input" << std::endl;
      return false;
    }
    if ( not lattice.Load(*child) )
    {
      std::cerr << "Error while reading Lattice from input" << std::endl;
      return false;
    }
 
    // loads lamarck functional
    child = handle.FirstChild( "Job" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find Functional in input" << std::endl;
      return false;
    }
    VA_CE::Functional_Builder ce;
    if ( not ce.Load(*child) )
    {
      std::cerr << "Error while reading Functional from input" << std::endl;
      return false;
    }
    ce.add_equivalent_clusters();
 
    // do structure
    child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
    for (; child; child = child->NextSiblingElement("Structure") )
    {
      Ising_CE::Structure structure;
      Ising_CE :: Structure :: lattice = &lattice;
      if ( not structure.Load(*child) )
      {
        std::cerr << "Error while reading Structure from input" << std::endl;
        return false;
      }
 
      VA_CE::Functional_Builder::t_VA_Functional functional;
      ce.generate_functional(structure, &functional);
      __MPICODE( functional.get_functional1()->set_mpi( &::mpi::main ); )
      __MPICODE( functional.get_functional2()->set_mpi( &::mpi::main ); )
    
      functional.resize( structure.atoms.size() );
      std::transform( structure.atoms.begin(), structure.atoms.end(), functional.begin(),
                      boost::lambda::bind( &Ising_CE::Structure::t_Atom::type,
                                           boost::lambda::_1 ) );
    
      std::cout << "Energy: " << functional.evaluate() << "\n"
                << "Concentration: " << structure.get_concentration() << "\n\n";
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
    __MPICODE( sstr << "\nProcessor " << mpi::main.rank() + 1
                    << " of " << mpi::main.size()
                    << " says:\n"; )

    sstr << "Caught error while running " << __PROGNAME__ 
         << "\n" << e.what() << "\n";

    __MPICODE( 
      mpi::AllGather bc( mpi::main );
      std::string message = sstr.str();
      bc << mpi::BroadCast::clear
         << message
         << mpi::BroadCast::allocate
         << message
         << mpi::BroadCast::broadcast;
      sstr.str("");
      while( bc.serialize( message ) ) sstr << message;
      bc << mpi::BroadCast::clear;
      ::mpi::main.barrier();
    )
    __NOTMPIROOT( return 0; )
    std::cerr << sstr.str() << std::endl;
    return 0;
  }
  return 1;
}

