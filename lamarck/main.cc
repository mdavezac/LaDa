//
//  Version: $Id$
//

#include <fstream>
#include <sstream>
#include <string>

#include <print/manip.h>

#include "functional_builder.h"

int main(int argc, char *argv[]) 
{
  TiXmlElement *child;
  atat::rVector3d vec;
  Ising_CE::Lattice lattice;

#ifdef _MPI
  mpi::main(argc, argv);
#endif
  
  
  std::string filename("input.xml");
  if( argc > 1 )
  {
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
               and (is_op == "-i" or is_op == "--input") ) istr >> filename;
      else if( is_op == "-h" or is_op == "--help" )
        std::cout << "Command-line options:\n\t -h, --help this message"
                  << "\n\t -i, --input XML input file (default: input.xml)\n\n";
    }
    filename = Print::reformat_home( filename );
    if( filename != "input.xml" )
      std::cout << "Reading from input file " << filename << std::endl;
  }

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
    std::cout << std::endl << std::endl;
    Ising_CE::Structure structure;
    Ising_CE :: Structure :: lattice = &lattice;
    if ( not structure.Load(*child) )
    {
      std::cerr << "Error while reading Structure from input" << std::endl;
      return false;
    }

    VA_CE::Functional_Builder::t_VA_Functional functional;
    ce.generate_functional(structure, &functional);
  
    functional.resize( structure.atoms.size() );
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = structure.atoms.end();
    VA_CE::Functional_Builder::t_VA_Functional::iterator i_var = functional.begin();
    for( ; i_atom != i_atom_end; ++i_atom, ++i_var )
      *i_var = i_atom->type;
  
    std::cout << "Energy: " << functional.evaluate() << std::endl
              << "Concentration: " << structure.get_concentration() << std::endl;
//   Ising_CE::Fourier( structure.atoms.begin(), structure.atoms.end(),
//                      structure.k_vecs.begin(), structure.k_vecs.end() );
//   structure.print_out( std::cout );

//   delete functional.get_functional1();
//   delete functional.get_functional2();
  }

}

