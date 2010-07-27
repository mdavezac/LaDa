#include "LaDaConfig.h"

#include <fstream>
#include <sstream>
#include <string>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <opt/mpi.h>
#include <opt/bpo_macros.h>

#include "functional_builder.h"

#define __PROGNAME__ "lamarck"

#define OUTPUT std::cout
#define ENDLINE "\n"
#ifndef TIXML_USE_STL
#error not using TIXML_USE_STL
#endif

#ifdef LADA_MPI
#include <boost/mpi/environment.hpp>
#endif

#include "constituent_strain.h"
#include "harmonic.h"

#if defined(_CUBIC_CE_)
typedef LaDa::CE::ConstituentStrain::Harmonic::Cubic t_Harmonic;
#elif defined( _TETRAGONAL_CE_ )
typedef LaDa::CE::ConstituentStrain::Harmonic::Tetragonal t_Harmonic;
#else
#error Please specify _CUBIC_CE_ or _TETRAGONAL_CE_
#endif
typedef LaDa::CE::Builder<t_Harmonic> t_Builder;

int main(int argc, char *argv[]) 
{
  namespace fs = boost::filesystem;
  LADA_TRY_BEGIN
  LADA_MPI_START

  __BPO_START__
  __BPO_RERUNS__;
  __BPO_HIDDEN__;
  po::options_description all; 
  all.add(generic);
  po::options_description allnhidden;
  allnhidden.add(all).add(hidden);
  po::positional_options_description p; 
  p.add("input", 1); 
  __BPO_MAP__
  __BPO_HELP_N_VERSION__;

  fs::path filename( vm["input"].as< std::string >() );
  LADA_DO_NASSERT(    ( not fs::exists( filename ) )
              or ( not ( fs::is_symlink(filename) or fs::is_regular(filename) ) ),
              filename << " does not exist.\n" )


  TiXmlElement *child;
  LaDa::math::rVector3d vec;
  LaDa::Crystal::Lattice lattice;
 
  
  LADA_MPI_CODE( std::string input; )
  LADA_ROOT( (*world),
    TiXmlDocument doc( filename.string() );
    LADA_DO_NASSERT( not doc.LoadFile(), 
                  "error while opening input file "
               << filename << "\n" << doc.ErrorDesc()  )
    LADA_MPI_CODE( 
      std::ostringstream stream;
      stream << *doc.RootElement();
      input = stream.str();
    )
  )
  
  LADA_MPI_CODE(
      boost::mpi::broadcast( *world, input, 0 ); 
      TiXmlDocument doc;
      doc.Parse( input.c_str() );
  )
  TiXmlHandle handle( &doc );
 
  // loads lattice
  child = handle.FirstChild( "Job" ).FirstChild( "Lattice" ).Element();
  LADA_DO_NASSERT( not child, "Could not find Lattice in input." )
  LADA_DO_NASSERT( not lattice.Load(*child), "Error while reading Lattice from input.")
 
  // loads lamarck functional
  child = handle.FirstChild( "Job" ).Element();
  LADA_DO_NASSERT( not child, "Could not find Functional in input." )

  t_Builder ce;
  LADA_DO_NASSERT( not ce.Load(*child), "Error while reading Functional from input." )
  ce.add_equivalent_clusters();
 
  LADA_ROOT( (*world), OUTPUT << "Nb procs: " << world->size() << ENDLINE; )
  // do structure
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  for (; child; child = child->NextSiblingElement("Structure") )
  {
    LaDa::Crystal::Structure structure;
    LaDa::Crystal :: Structure :: lattice = &lattice;
    LADA_DO_NASSERT( not structure.Load(*child), "Error while reading Structure from input." )
 
    t_Builder::t_VA_Functional functional;
    ce.generate_functional(structure, &functional);
    LADA_MPI_CODE( functional.get_functional1()->set_mpi( world.get() ); )
    LADA_MPI_CODE( functional.get_functional2()->set_mpi( world.get() ); )
  
    functional.resize( structure.atoms.size() );
    std::transform( structure.atoms.begin(), structure.atoms.end(), functional.begin(),
                    boost::lambda::bind( &LaDa::Crystal::Structure::t_Atom::type,
                                         boost::lambda::_1 ) );
  
    OUTPUT << "Energy: " << functional.evaluate() // << "\n";
           << "  Concentration: " << structure.get_concentration() << "\n" ; //<< ENDLINE;
// LaDa::Crystal::Fourier( structure.atoms.begin(), structure.atoms.end(),
//                    structure.k_vecs.begin(), structure.k_vecs.end() );
// structure.print_out( std::cout );

// delete functional.get_functional1();
// delete functional.get_functional2();
  }
  return 0;
  __BPO_CATCH__()
}

