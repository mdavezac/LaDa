#include "functional.h"
#include <lamarck/lattice.h>
#include <lamarck/structure.h>
#include <tinyxml/tinyxml.h>
#include <opt/opt_frprmn.h>
#include <opt/opt_minimize_gsl.h>


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
  for (; child ; child = child->NextSiblingElement("Structure") )
  {
    if ( not structure.Load(*child) )
    {
      std::cerr << "Error while reading Lattice from input" << std::endl;
      return false;
    }

    // loads vff functional
    TiXmlElement *vff_xml = handle.FirstChild( "Job" ).Element();
    if ( not vff_xml )
    {
      std::cerr << "Could not find Lattice in input" << std::endl;
      return false;
    }
    Vff::Functional vff(structure);
    if ( not vff.Load(*vff_xml) )
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
    std::cout << "Energy: " << structure.energy << std::endl;
    std::cout << "Stress Tensor: " << std::endl;
    const atat::rMatrix3d stress = vff.get_stress();
    std::cout << "   " << stress(0,0) << " " << stress(1,0) << " " << stress(2,0) << std::endl
              << "   " << stress(0,1) << " " << stress(1,1) << " " << stress(2,1) << std::endl
              << "   " << stress(0,2) << " " << stress(1,2) << " " << stress(2,2) 
              << std::endl << std::endl << std::endl;

  }

  return 0;
}
