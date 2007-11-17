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
  mpi::main(argc, argv);

  std::ostringstream sstr; sstr << "Am Proc " << mpi::main.rank();
  std::string str = sstr.str();
  if( not mpi::main.is_root_node() ) std::cout << str << std::endl;
  mpi::BroadCast bc(mpi::main);
  bc << str << mpi::BroadCast::allocate << str << mpi::BroadCast::broadcast << str;
  if( not mpi::main.is_root_node() ) std::cout << str << std::endl;
  return 1;
}
