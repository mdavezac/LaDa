//
//  Version: $Id$
//
// #include <tinyxml/tinyxml.h>
//

#include <sstream>
#include <string>

#include "mpi_object.h"


int main(int argc, char *argv[]) 
{
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  ::mpi::main = &world;

//  std::cout << "rank: " << mpi::main.rank() << std::endl;
//  std::ostringstream sstr; sstr << "Am Proc " << mpi::main.rank();
//  std::string str = sstr.str();
//  if( not ::mpi::main.is_root_node() ) str = "";
//  mpi::BroadCast bc(mpi::main);
//  bc.serialize( str );
//  bc.allocate_buffers();
//  bc.serialize( str );
//  bc();
//  bc.serialize( str );
// / bc << str << mpi::BroadCast::allocate << str << mpi::BroadCast::broadcast << str;
//  mpi::main.barrier();
//  if( not ::mpi::main.is_root_node() ) std::cout << "end" << std::endl;
//  mpi::main.barrier();
//  if( ::mpi::main.rank() ) std::cout << std::endl << " done \n" << str << std::endl; 
  return 1;
}
