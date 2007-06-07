
#include <fstream>
#include <sstream>
#include <string>

#include <lamarck/functional_builder.h>

int main(int argc, char *argv[]) 
{
  TiXmlElement *child;
  atat::rVector3d vec;
  Ising_CE::Lattice lattice;
  
  
  std::string filename("input.xml");
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
    structure.set_kvectors( functional.get_functional2()->get_kvectors() );
  
    functional.resize( structure.atoms.size() );
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = structure.atoms.end();
    VA_CE::Functional_Builder::t_VA_Functional::iterator i_var = functional.begin();
    for( ; i_atom != i_atom_end; ++i_atom, ++i_var )
      *i_var = i_atom->type;
  
    std::cout << "Energy: " << functional.evaluate() << std::endl
              << "Concentration: " << structure.get_concentration() << std::endl;
    Ising_CE::fourrier_to_kspace( structure.atoms.begin(), structure.atoms.end(),
                                  structure.k_vecs.begin(), structure.k_vecs.end() );
    structure.print_out( std::cout );
  }

}

