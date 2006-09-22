#include "motu.h"
using LaDa::MotU;


int main(int argc, char *argv[]) 
{
  MotU motu("input.xml");
  motu.run();
  return 0;
}
