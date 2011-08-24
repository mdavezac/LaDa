#include "LaDaConfig.h"
#include<iostream>
#include<string>
#include<vector>

#include "../supercell.h"
#include "../primitive.h"

using namespace std;
int main()
{
  using namespace LaDa;
  using namespace LaDa::crystal;
  using namespace LaDa::math;
  TemplateStructure< std::vector<std::string> > lattice; 
  lattice.set_cell(0,0.5,0.5)
                  (0.5,0,0.5)
                  (0.5,0.5,0);
  lattice.add_atom(0,0,0, "Si")
         .add_atom(0.5,0.5,0.5, "Si", "Ge");
  rMatrix3d cell;
  size_t nmax= 8;
  bool cont = true;
  for(size_t a(1); a <= nmax and cont; ++a)
  {
    if(nmax % a != 0) continue;
    size_t const Ndiv_a = nmax/a;
    cell(0,0) = a;
    // Iterates over values of a such that a * b * c == nmax
    for(size_t b(1); b <= Ndiv_a and cont; ++b) 
    {
      cell(1,1) = b;
      if(Ndiv_a % b != 0) continue;
      // Iterates over values of a such that a * b * c == nmax
      size_t const c( Ndiv_a/b);
      cell(2,2) = c;
      if( a * b *c != nmax ) BOOST_THROW_EXCEPTION(LaDa::error::internal())
      for(size_t d(0); d < b and cont; ++d) 
      {
        cell(1,0) = d;
        for(size_t e(0); e < c and cont; ++e) 
        {
          cell(2,0) = e;
          for(size_t f(0); f < c and cont; ++f) 
          {
            cell(2,1) = f;
            structure = supercell(lattice, cell);
            std::cout << structure << "\n";
            LADA_DOASSERT(is_integer(structure.cell() * lattice.cell().inverse()), "Not a sublattice.");
            LADA_DOASSERT(is_integer(lattice.cell() * structure.cell().inverse()), "Not a sublattice.");
            LADA_DOASSERT(eq(lattice.cell().determinant(), structure.cell().derterminant()), "Not primitive.");
            cont = false;
          } // f
        } // e
      } // d
    } // b
  } // a

  return 0;
}
