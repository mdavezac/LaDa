//
//  Version: $Id$
//

#include <tinyxml/tinyxml.h>

#include "functional.h"
#include <lamarck/lattice.h>
#include <lamarck/structure.h>
#ifdef _DOFORTRAN
#include <opt/opt_frprmn.h>
#else
#include <opt/gsl_minimizers.h>
#endif


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
    
#ifdef _DOFORTRAN
    Minimizer::Frpr<Vff::Functional> minimizer( vff, vff_for_frprmn );
    vff_for_frprmn( (double*) &vff, NULL );
#else
    Minimizer::GnuSL<Vff::Functional> minimizer( vff );
#endif
    child = handle.FirstChild( "Job" ).Element();
    minimizer.Load(*child);
    minimizer();
    structure.energy = vff.energy() / 16.0217733;
    const atat::rMatrix3d stress = vff.get_stress();
    std::cout << std::fixed << std::setprecision(5) 
              << "Energy [eV]: " << std::setw(12) << structure.energy << std::endl
              << "Stress Tensor: " << std::endl 
              << std::setw(12) << stress(0,0) << " " << std::setw(12)
                               << stress(1,0) << " " << std::setw(12) << stress(2,0) << std::endl
              << std::setw(12) << stress(0,1) << " " << std::setw(12)
                               << stress(1,1) << " " << std::setw(12) << stress(2,1) << std::endl
              << std::setw(12) << stress(0,2) << " " << std::setw(12)
                               << stress(1,2) << " " << std::setw(12) << stress(2,2) 
              << std::endl << std::endl << std::endl;
    vff.print_escan_input( "atomic.config" );


    Ising_CE::Fourier( structure.atoms.begin(),  structure.atoms.end(),
                       structure.k_vecs.begin(), structure.k_vecs.end() );
    structure.print_out(std::cout );

  }

  return 0;
}
