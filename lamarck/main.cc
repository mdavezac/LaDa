//
//  Version: $Id$
//

#include <fstream>
#include <sstream>
#include <string>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <print/manip.h>

#include "functional_builder.h"

#include <revision.h>
#define __PROGNAME__ "lamarck"

void parse_cli( int argc, char *argv[], std::string &_filename )
{
  if( argc < 2 ) return;
  
  std::ostringstream sstr;
  for( types::t_int i = 1; i < argc; ++i )
    sstr << argv[i] << " "; 
  std::istringstream istr( sstr.str() );
  while ( istr.good() )
  {
    std::string is_op;
    istr >> is_op; is_op = Print::StripEdges( is_op );
    if( is_op.empty() ) continue;
    else if(     istr.good()
             and (is_op == "-i" or is_op == "--input") ) istr >> _filename;
    else if( is_op == "-h" or is_op == "--help" )
    {
      std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package.\n"
                << "Command-line options:\n\t -h, --help this message"
                << "\n\t -v, --version Subversion Revision and Package version"
                << "\n\t -i, --input XML input file (default: input.xml)\n\n";
      exit(1);
    }
    else if( is_op == "-v" or is_op == "--version" )
    {
      std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package\n"
                << "Subversion Revision: " << SVN::Revision << "\n\n"; 
      exit(1);
    }
  }
  _filename = Print::reformat_home( _filename );
}

int main(int argc, char *argv[]) 
{
  TiXmlElement *child;
  atat::rVector3d vec;
  Ising_CE::Lattice lattice;

  std::string filename("input.xml");
  __MPICODE( mpi::main(argc, argv); )
  __ROOTCODE( parse_cli( argc, argv, filename ); )
  
  
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

