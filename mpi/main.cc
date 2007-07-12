#include <tinyxml/tinyxml.h>

#include "lamarck/functional_builder.h"
#include "lamarck/structure.h"
#include "vff/functional.h"
#include "opt/opt_minimize_gsl.h"
#include "mpi_object.h"


int main(int argc, char *argv[]) 
{
  TiXmlElement *child;
  atat::rVector3d vec;
  VA_CE :: Functional_Builder builder;
  VA_CE :: Functional_Builder :: t_VA_Functional functional;
  Ising_CE::Structure structure;
  
  mpi::main(argc, argv);

  std::string filename("input.xml");
  TiXmlDocument doc( filename.c_str() );
  TiXmlHandle handle( &doc );
  if ( mpi::main.rank() == 0 )
  {
    if  ( !doc.LoadFile() )
    {
      std::cerr << "error while opening input file " << filename << std::endl
                << doc.ErrorDesc() << std::endl; 
      return false;
    }
    
    const TiXmlElement *functional_xml = handle.FirstChild("Job")
                                               .FirstChild("Functional").Element();
    for(; functional_xml; functional_xml = functional_xml->NextSiblingElement("Functional") )
    {
      std::string str = ""; 
      if ( functional_xml->Attribute("type") )
        str = functional_xml->Attribute("type");
      if ( str.compare("CE") == 0 )
        break;
    }
    if ( not functional_xml )
    {
      std::cerr << "No <Functional type=\"CE\"> found "
                << std::endl << "Giving up" << std::endl;
      return false;
    }
 
    if ( not builder.Load(*functional_xml) )
      return false;
 
    // loads structure
    child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find Lattice in input" << std::endl;
      return false;
    }
    if ( not structure.Load(*child) )
    {
      std::cerr << "Error while reading Lattice from input" << std::endl;
      return false;
    }
  }
 
  types::t_real *buff = new types::t_real[10*mpi::main.size()];
  types::t_real *uff = buff;
  for( types::t_int i = 0; i < 10; ++i, ++uff)
    *uff = i + ( pow(10,mpi::main.rank()) );
 
  mpi::AllGather allgather( mpi::main );
  allgather.serialize( buff, buff+10);
  allgather.allocate_buffers(0);
  allgather.serialize( buff, buff+10);
  allgather();
  allgather.serialize( buff, buff+10*mpi::main.size());
 
  for( types::t_int i = 0; i < 10 * mpi::main.size(); ++i)
    std::cout << *(buff+i) << std::endl;
  mpi::BroadCast broadcast( mpi::main );
  broadcast.serialize(structure);
  broadcast.serialize(builder);
  broadcast.allocate_buffers();
  broadcast.serialize(structure);
  broadcast.serialize(builder);
  broadcast();
  broadcast.serialize(structure);
  broadcast.serialize(builder);
 
  builder.add_equivalent_clusters();
 
  builder.generate_functional(structure, &functional);
  structure.set_kvectors( functional.get_functional2()->get_kvectors() );
  functional.resize( structure.atoms.size() );
  VA_CE::Functional_Builder::t_VA_Functional::iterator i_var = functional.begin();
  VA_CE::Functional_Builder::t_VA_Functional::iterator i_var_end = functional.begin();
  Ising_CE::Structure::t_Atoms::iterator i_atom = structure.atoms.begin();
  for(; i_var != i_var_end; ++i_var, ++i_atom )
    *i_var = i_atom->type;
 
  if ( mpi::main.rank() == 0 )
    std::cout << " Size: " << mpi::main.size() << std::endl << std::endl;
  for( int i=0; i < mpi::main.size(); ++i)
  {
    if ( i == mpi::main.rank() )
    {
      std::cout << " Rank: " <<  i << " " << functional.evaluate() << std::endl;
      structure.print_out( std::cout );
      functional.get_functional1()->print_out(std::cout);
      std::cout << std::endl;
    }
    mpi::main.barrier();
  }

  delete functional.get_functional1();
  delete functional.get_functional2();

  return 0;
}
