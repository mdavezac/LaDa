//
//  Version: $Id$
//
#ifndef _VFF_FUNCTIONAL_H_
#define _VFF_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <algorithm>
#include <boost/filesystem/path.hpp>

#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <opt/types.h>
#include <opt/function_base.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

#include "atomic_center.h"
#include "atomic_functional.h"

#ifdef _MPI
# include <mpi/mpi_object.h>
#endif 

namespace LaDa
{
  //! \cond
  namespace Crystal { template< class T_TYPE > class ConquerBox; }
  //! \endcond 


  //! \brief Reimplements the Valence Force Field %Functional in c++
  //! \details Vff, or Valence Force Field functional, is an empirical functional which
  //! attempts to model the strain energy of a material from at most three-body
  //! interactions. The there body interactions which are considered are
  //! bond-stretching, change in bond-angles, and a combination of these two.
  //! 
  //! The implementation relies on a body-centered paradigm. In other words,
  //! four classes have been created:
  //!   - Vff::Vff is wrapper class and interface to Vff
  //!   - Vff::AtomicCenter represent a single atom and lists its first neighbor relationships
  //!   - Vff::AtomicCenter::const_iterator allows coders to travel
  //!   through a collection of Vff::Atomic_Centera along first neighbor
  //!   relationships.
  //!   - Vff::AtomicFunctional computes the strain energy of one single atom, eg
  //!   all three body terms in which a particular atom takes part.
  //!   .
  namespace Vff
  {
    //! \brief Wrapper, or interface class to Valence Force Field Vff
    //! \details Takes care of loading parameters from input, building tree for Functional::structure,
    //! and eventually is able to compute its the energy and strain.
    //! One possible use is the following:
    //! \code
    //! Crystal::Structure structure;
    //! Crystal::Lattice lattice;
    //! Crystal :: Structure :: lattice = &lattice; // don't forget this hack
    //! // load lattice, structure, etc from input
    //! // Then load Vff input itself
    //! TiXmlElement *vff_xml = handle.FirstChild( "Job" ).Element();
    //! Vff::Functional vff(structure);
    //! if ( not vff.Load(*vff_xml) ) return false;
    //! 
    //! vff.initialize_centers(); // Construct mesh of first neighbor relationships
    //!
    //! 
    //! Minimizer::GnuSL<Vff::Functional> minimizer( vff ); // Declare and Load a GSL minimizer
    //! child = handle.FirstChild( "Job" ).Element();
    //! minimizer.Load(*child);
    //!
    //! minimizer.minimize(); // Minimize strain
    //! structure.energy = vff.energy() / 16.0217733; // compute relaxed strain energy
    //! std::cout << std::fixed << std::setprecision(5)  // prints stuff out
    //!           << "Energy [meV/atom]: " << std::setw(12) << structure.energy << std::endl
    //!           << "Stress Tensor: " << std::endl 
    //!           << std::setw(12) << stress(0,0) << " " << std::setw(12)
    //!                            << stress(1,0) << " " << std::setw(12) << stress(2,0) << std::endl
    //!           << std::setw(12) << stress(0,1) << " " << std::setw(12)
    //!                            << stress(1,1) << " " << std::setw(12) << stress(2,1) << std::endl
    //!           << std::setw(12) << stress(0,2) << " " << std::setw(12)
    //!                            << stress(1,2) << " " << std::setw(12) << stress(2,2) 
    //!           << std::endl << std::endl << std::endl;
    //! vff.print_escan_input( "atomic.config" ); // print pescan input
    //! 
    //! \endcode
    class Vff __MPICODE( : public MPI_COMMDEC )
    {
      protected:
        //! Type of the path.
        typedef boost::filesystem::path t_Path;
        //! The type of the atom  
        typedef Crystal::Structure::t_Atom  t_Atom;
        //! The type of the atom container
        typedef Crystal::Structure::t_Atoms t_Atoms;

      public:
        //! \brief Constructor and Initializer
        //! \param _str structure for which to compute energy and stress
        Vff   ( Crystal :: Structure &_str )
            : structure(_str),
              bond_cutoff(0) {}
        //! \brief Copy Constructor
        Vff   ( const Vff &_c )
            : __MPICODE( MPI_COMMCOPY( _c ) __COMMA__ )
              structure( _c.structure ),
              bond_cutoff( _c.bond_cutoff ),
              centers( _c.centers ), functionals( _c.functionals ) {}
        //! \brief Destructor
        ~Vff() {}

        //! \brief Loads input to functional from  xml 
        //! \param _element should point to an xml node which is the functional data
        //! or contains a node to the funtional data as a direct child
        bool Load( const TiXmlElement &_element );

        //! \brief computes energy and stress, expects everything to be set
        types::t_real energy() const;
        //! \brief Prints atom.config type input to escan
        //! \param _f optional filename to which to direct the output
        void print_escan_input( const t_Path &_f = "atom.config") const;
        //! \deprecated Constructs the mesh of AtomicCenter
        //! Defined in initialize_centers.cc.
        bool initialize_centers();
        //! Prints out all parameters
        void print_out( std::ostream &stream ) const;

      protected:
        //! Holds ideal first neighbor positions.
        typedef std::vector< std::vector< atat::rVector3d > > t_FirstNeighbors;
        //! Type of the smith normal transformation.
        typedef boost::tuples::tuple
                < 
                  atat::rMatrix3d, 
                  atat::iVector3d 
                > t_Transformation;
        
        //! \brief Loads Functional for one or two site lattices.
        //! \details If \a _node is not the correct node, the results are undefined.
        bool load_( const TiXmlElement &_node );
        //! Finds first neighbors of ideal lattice sites.
        void first_neighbors_( t_FirstNeighbors& _fn );

        //! \brief Constructs the mesh of AtomicCenter
        //! \details Newer version than Functional::construct_centers. 
        //!          Uses smith normal form to speed up first neighbor search.
        //!          Defined in build_tree.cc.
        bool build_tree_smith_(const t_FirstNeighbors& _fn);
        //! Builds first neighbor tree using two loops over centers.
        bool build_tree_sort_(const t_FirstNeighbors& _fn);
        //! Builds first neighbor tree using two loops over centers and divide-n-conquer.
        bool build_tree_sort_dnc_( const Crystal::ConquerBox<types::t_real>& _dnc, 
                                   const t_FirstNeighbors& _fn);

        //! \brief computes smith index.
        //! \details Defined in build_tree.cc.
        void smith_index_( const t_Transformation &_transformation,
                           const atat::rVector3d &_pos,
                           atat::iVector3d &_index );
        //! \brief Computes matrix to get to smith index.
        //! \details Defined in build_tree.cc.
        t_Transformation to_smith_matrix( const atat::rMatrix3d &_lat_cell,
                                          const atat::rMatrix3d &_str_cell );

        //! Type of the atomic centers
        typedef AtomicCenter t_Center;  
        //! Type of the container holding the atomic centers
        typedef t_Center :: t_Centers t_Centers;  
        //! Type of the container holding the atomic functionals
        typedef std::vector< AtomicFunctional > t_AtomicFunctionals;  
        
        //! Crystal::Structure for which to compute energy and stress
        mutable Crystal :: Structure &structure;
        //! length below which first-neighbor relationship is defined
        types::t_real bond_cutoff; 
        //! \brief list of all AtomicCenter created from Vff::structure
        //! \details Space for the atomic centers are reserved in the
        //!          constructor. It is expected that the iterators will be valid
        //!          throughout the life of the functional, starting with a call
        //!          to the construction of the tree. If the order of the atoms
        //!          in Vff::structure is changed, then it is imperative
        //!          that the tree be reconstructed from scratch.
        t_Centers centers;  
        //! list of all possbile Atomic_Functionals for Vff::structure.lattice
        t_AtomicFunctionals functionals;
        

#      ifdef _LADADEBUG
         //! Checks that the list of centers are valid. somewhat.
         void check_tree() const;
#      endif
    };

  } // namespace vff 
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
