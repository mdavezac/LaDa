#ifndef _VFF_FUNCTIONAL_H_
#define _VFF_FUNCTIONAL_H_

#include "LaDaConfig.h"

#include <vector>
#include <algorithm>
#include <boost/filesystem/path.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/int.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/dynamic_bitset.hpp>


#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <crystal/periodic_dnc.h>
#include <misc/types.h>
#include <opt/function_base.h>
#include <opt/debug.h>
#include <opt/mpi.h>
#include <math/eigen.h>

#include "atomic_center.h"
#include "data.h"
#include "exceptions.h"

#ifdef LADA_MPI
# include <opt/mpi.h>
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
  namespace vff
  {
    //! Valence Force Field functional. 
    class Vff 
    {
      public:
        //! Return of Vff.
        typedef types::t_real t_Return;
        //! Argument to Vff.
        typedef Crystal::TStructure<std::string> t_Arg;
        //! Keeps track of structure being currently evaluated.
        mutable t_Arg structure;
        //! First-neighbor bond parameter.
        types::t_real bond_cutoff;
        //! Bond parameters.
        std::map<std::string, BondData> bonds_params;
        //! Angle parameters.
        std::map<std::string, AngleData> angles_params;
#       ifdef LADA_MPI
          //! MPI communicator
          boost::mpi::communicator comm;
#       endif
        
      protected:
        //! Type of the path.
        typedef boost::filesystem::path t_Path;
        //! The type of the atom  
        typedef t_Arg::t_Atom  t_Atom;
        //! The type of the atom container
        typedef t_Arg::t_Atoms t_Atoms;

      public:
        //! Constructor and Initializer
        Vff() : bond_cutoff(0) {}
        //! \brief Constructor and Initializer
        //! \param _str structure for which to compute energy and stress
        Vff(t_Arg const &_str) : structure(_str), bond_cutoff(0) {}
        //! \brief Copy Constructor
        Vff   ( const Vff &_c )
            : structure( _c.structure ), bond_cutoff( _c.bond_cutoff ), 
              centers_( _c.centers_ ), angles_params(_c.angles_params), bonds_params(_c.bonds_params) {}
        //! \brief Destructor
        ~Vff() {}

#       ifdef LADA_MPI
          //! \brief computes energy and stress, expects everything to be set
          t_Return operator()(t_Arg const &_arg, boost::mpi::communicator const &_comm)
            { comm = _comm; return operator()(_arg); }
          //! \brief computes energy and stress, expects everything to be set
          t_Return energy(boost::mpi::communicator const &_comm) 
            { comm = _comm; return energy(); }
#       endif
        //! \brief computes energy and stress, expects everything to be set
        t_Return operator()(t_Arg const &_arg)
        { 
          structure = _arg;
          initialize_centers();
          return energy();
        }
        //! \brief computes energy and stress, expects everything to be set
        types::t_real energy() const;
        //! \brief Prints atom.config type input to escan
        //! \param _f optional filename to which to direct the output
        void print_escan_input( const t_Path &_f = "atom.config") const;
        //! \brief Constructs the mesh of AtomicCenter
        //! \details Defined in initialize_centers.cc.
        bool initialize_centers(bool _verbose = false);

        //! \brief checks that all required bond/angle parameters are there. 
        //! \details throws if input is missing.
        void check_input() const;

      protected:
        //! Holds ideal first neighbor positions.
        typedef std::vector< std::vector< math::rVector3d > > t_FirstNeighbors;
        //! Type of the smith normal transformation.
        typedef boost::tuples::tuple
                < 
                  math::rMatrix3d, 
                  math::iVector3d 
                > t_Transformation;
        
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
        bool build_tree_sort_dnc_( const Crystal::ConquerBox<t_Atom::t_Type>& _dnc, 
                                   const t_FirstNeighbors& _fn);
        //! Builds first neighbor tree by sorting distances in divide and conquer box.
        bool build_tree_partial_sort_dnc( const Crystal::DnCBoxes::value_type& _dnc, 
                                          types::t_real _cutoff);

        //! \brief computes smith index.
        //! \details Defined in build_tree.cc.
        void smith_index_( const t_Transformation &_transformation,
                           const math::rVector3d &_pos,
                           math::iVector3d &_index );
        //! \brief Computes matrix to get to smith index.
        //! \details Defined in build_tree.cc.
        t_Transformation to_smith_matrix( const math::rMatrix3d &_lat_cell,
                                          const math::rMatrix3d &_str_cell );

        //! Returns bond data for given bond.
        BondData const & get_bond_data(AtomicCenter::const_iterator const _i_bond) const
          { return bonds_params.find(_i_bond.key())->second; }
        //! Returns angle data for given angle.
        AngleData const & get_angle_data( AtomicCenter::const_iterator const _i_A,
                                          AtomicCenter::const_iterator const _i_B ) const
          { return angles_params.find(_i_A.angle_key(_i_B))->second; }

        //! Computes micro strain around given atom.
        types::t_real micro_strain(const AtomicCenter &_center) const;
        //! Evaluates VFF around single atom.
        types::t_real evaluate_center(const AtomicCenter &_center) const;
        //! \brief Evaluate strain energy and gradients for AtomicCenter _center
        //! \details returns strain energy, and computes stress
        //! \param _center center for which to evaluate energy
        //! \param _strain to which AtomicFunctional::structure is submitted
        //! \param _stress on _center resulting from _strain
        //! \param _K0 displacement resulting from _strain and w.r.t original unit-cell
        //! \sa function::Base, function::Base::evaluate_with_gradient()
        types::t_real evaluate_center_with_gradient( const AtomicCenter &_center,
                                                     const math::rMatrix3d &_strain,
                                                     math::rMatrix3d &_stress,
                                                     const math::rMatrix3d &_K0 ) const;
        //! Type of the atomic centers
        typedef AtomicCenter t_Center;  
        //! Type of the container holding the atomic centers
        typedef t_Center :: t_Centers t_Centers;  
        
        //! \brief list of all AtomicCenter created from Vff::structure
        //! \details Space for the atomic centers are reserved in the
        //!          constructor. It is expected that the iterators will be valid
        //!          throughout the life of the functional, starting with a call
        //!          to the construction of the tree. If the order of the atoms
        //!          in Vff::structure is changed, then it is imperative
        //!          that the tree be reconstructed from scratch.
        t_Centers centers_;  
    };
  } // namespace vff 
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
