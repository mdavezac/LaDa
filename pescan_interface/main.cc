#include <vff/functional.h>
#include <lamarck/lattice.h>
#include <lamarck/structure.h>
#include <tinyxml/tinyxml.h>
#include <opt/opt_frprmn.h>
#include <opt/opt_minimize_gsl.h>

#include "interface.h"

extern "C" {
  double call_it( const double* const _x,
                  double* const _y)
  {
    static Vff::Functional* this_func;
    if ( not _y  )
    {
      this_func = (Vff::Functional*) _x;
      return 0;
    }

    if ( not this_func )
      return 0;

    const double *i_x_copy = _x;
    const double *i_x_end = _x + this_func->size();
    std::copy(i_x_copy, i_x_end, this_func->begin());
    types::t_real result = this_func->evaluate_with_gradient(_y);
    std::cout << "Energy " << result << std::endl;
    return result;
  }
}

int main(int argc, char *argv[]) 
{
#ifdef _MPI
  mpi::main(argc, argv);
#endif
  TiXmlElement *child;
  atat::rVector3d vec;
  Ising_CE::Structure structure;
  Ising_CE::Lattice lattice;
  
  Ising_CE :: Structure :: lattice = &lattice;
  
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

  // loads structure
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  if ( not child )
  {
    std::cerr << "Could not find Lattice in input" << std::endl;
    return false;
  }
  if ( not structure.Load(*child) )
  {
    std::cerr << "Error while reading Structure from input" << std::endl;
    return false;
  }

  // loads vff functional
  child = handle.FirstChild( "Job" ).Element();
  if ( not child )
  {
    std::cerr << "Could not find Vff in input" << std::endl;
    return false;
  }
  Vff::Functional vff(structure);
  if ( not vff.Load(*child) )
  {
    std::cerr << "Error while reading Lattice from input" << std::endl;
    return false;
  }

  vff.initialize_centers();
  minimizer::GnuSL<Vff::Functional> minimizer( vff );
  child = handle.FirstChild( "Job" ).Element();
  minimizer.Load(*child);
  minimizer.minimize();
  structure.energy = vff.energy();
  vff.print_escan_input();

  Pescan::Interface interface;
  interface.Load(*handle.FirstChild("Job").Element());
  std::cout << "Band Gap: " << interface(structure) << std::endl;

  return 0;
}
