

#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>
#include "../structure.h"
#include "../is_container.h"

int main()
{
  LaDa::crystal::TemplateStructure<std::string> structure;
  structure.set_cell(-0.5,0.5,0.5)
                    (0.5,-0.5,0.5)
                    (0.5,0.5,-0.5);
  structure.add_atom(0,0,0,"Au");

  std::cout << structure << "\n";

  return 0;
};
