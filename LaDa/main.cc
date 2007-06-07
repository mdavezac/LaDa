#include "lamarck.h"
#include "darwin.h"

int main(int argc, char *argv[]) 
{
  try
  {
    LaDa::Lamarck lamarck;
    if ( not  lamarck.Load("input.xml") )
      throw ""; 
    LaDa::Darwin<LaDa::Lamarck::t_Individual, LaDa::Lamarck> darwin(  &lamarck );
    if ( not  darwin.Load("input.xml") )
      throw ""; 
    darwin.run();
  }
  catch ( std::exception &e )
  {
    std::cerr << "Caught error while running lada" << std::endl
	      << e.what();
  }
  return 0;
}
