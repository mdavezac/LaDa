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


#include <tinyxml/tinyxml.h>

#include <crystal/atom.h>
#include <crystal/structure.h>
#include <crystal/lattice.h>
#include <opt/types.h>
#include <opt/function_base.h>
#include <opt/debug.h>
#include <opt/mpi.h>

#include "atomic_center.h"
#include "atomic_functional.h"

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
    class Vff LADA_MPI_CODE( : public MPI_COMMDEC )
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
            : LADA_MPI_CODE( MPI_COMMCOPY( _c ) LADA_COMMA )
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
        //! \brief Constructs the mesh of AtomicCenter
        //! \details Defined in initialize_centers.cc.
        bool initialize_centers(bool _verbose = false);
        //! Prints out all parameters
        void print_out( std::ostream &stream ) const;

        //! Sets the bond parameters.
        template< class T_TUPLE >
          void set_bond( const std::string &_type, const T_TUPLE& _tuple );
        //! Returns bond parameters, first the length, then the alphas.
        boost::tuples::tuple< const types::t_real&, const types::t_real&,
                              const types::t_real&, const types::t_real&, 
                              const types::t_real&, const types::t_real& >
          get_bond( const std::string &_type ) const;
        //! Sets the angle parameters.
        template< class T_TUPLE >
          void set_angle( const std::string &_type, const T_TUPLE& _tuple );
        //! Returns angle parameters, first the length, then sigma, then the betas.
        boost::tuples::tuple< const types::t_real&, const types::t_real&, 
                              const types::t_real&, const types::t_real&,
                              const types::t_real&, const types::t_real&,
                              const types::t_real& >
          get_angle( const std::string &_type ) const;

        //! Copies parameters from argument.
        void copy_parameters(Vff const &_f)
          { functionals = _f.functionals; bond_cutoff = _f.bond_cutoff; }

      protected:
        //! Holds ideal first neighbor positions.
        typedef std::vector< std::vector< math::rVector3d > > t_FirstNeighbors;
        //! Type of the smith normal transformation.
        typedef boost::tuples::tuple
                < 
                  math::rMatrix3d, 
                  math::iVector3d 
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
                           const math::rVector3d &_pos,
                           math::iVector3d &_index );
        //! \brief Computes matrix to get to smith index.
        //! \details Defined in build_tree.cc.
        t_Transformation to_smith_matrix( const math::rMatrix3d &_lat_cell,
                                          const math::rMatrix3d &_str_cell );

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


#      ifdef LADA_DEBUG
         //! Checks that the list of centers are valid. somewhat.
         void check_tree() const;
#      endif
    };

    namespace details
    {
      struct type_index
      {
        int operator()( const Crystal::Lattice::t_Site::t_Type& _type,
                        const std::string & _s ) const
        {     
          typedef Crystal::Lattice::t_Site::t_Type t_Type;
          const t_Type :: const_iterator i_found
            = std::find( _type.begin(), _type.end(), _s );
          if( i_found == _type.end() ) return -1;
          return (int)( i_found - _type.begin() );
        }
      };
    }


    template< class T_TUPLE >
      void Vff::set_bond( const std::string &_type, const T_TUPLE& _tuple )
      {
        namespace bt = boost::tuples;

        LADA_DO_NASSERT( not (structure.lattice and structure.lattice->sites.size() == 2),
                    "Lattice undefined or does not have two sites.\n" )
        namespace bx = boost::xpressive;
        const std::string bond( boost::algorithm::trim_copy( _type ) );
        bx::smatch what;

        const bx::sregex regex =(    ( bx::s1 = ( bx::alpha >> !bx::alpha ) )
                                  >> *bx::_s >> "-" >> *bx::_s
                                  >> ( bx::s2 = ( bx::alpha >> !bx::alpha ) ) );
        if( !bx::regex_match( bond, what, regex ) )
          LADA_DO_NASSERT( true, "Could not find bond " + _type + "\n" )

        const bool two_species( structure.lattice->sites[0].type.size() == 2 );
        const std::string A( what.str(1) );
        const std::string B( what.str(2) );
        for( size_t site(0); site < 2; ++site )
        { 
          for( size_t swap(0); swap < 2; ++swap )
          {
            const int
              i = details::type_index()( structure.lattice->sites[ site ].type,  swap ? A: B ),
              j = details::type_index()( structure.lattice->sites[ site ? 0: 1 ].type,  swap ? B: A );
            if( i == -1 or j == -1 ) continue;
            const size_t Akind
            (
              site == 0 ? size_t(i):( two_species ? size_t(2 + i): size_t(1 + i) )
            );
            const size_t bond_kind
            (
              site == 0 ? size_t(j):( two_species ? size_t(j): 0 )
            );
            LADA_DO_NASSERT( Akind >= functionals.size(), "Index out-of-range.\n" )
            functionals[Akind].set_bond( bond_kind, _tuple );
          }
        }
      }
    
    template< class T_TUPLE >
      void Vff::set_angle( const std::string &_type, const T_TUPLE& _tuple )
      {
        namespace bt = boost::tuples;

        LADA_DO_NASSERT( not (structure.lattice and structure.lattice->sites.size() == 2),
                    "Lattice undefined or does not have two sites.\n" )
        namespace bx = boost::xpressive;
        const std::string bond( boost::algorithm::trim_copy( _type ) );
        bx::smatch what;

        const bx::sregex regex =(    ( bx::s1 = ( bx::alpha >> !bx::alpha ) )
                                  >> *bx::_s >> "-" >> *bx::_s
                                  >> ( bx::s2 = ( bx::alpha >> !bx::alpha ) )
                                  >> *bx::_s >> "-" >> *bx::_s
                                  >> ( bx::s3 = ( bx::alpha >> !bx::alpha ) ) );
        if( !bx::regex_match( bond, what, regex ) )
          LADA_DO_NASSERT( true, "Could not find angle " + _type + "\n" )

        const bool two_species( structure.lattice->sites[0].type.size() == 2 );
        const std::string A( what.str(1) );
        const std::string B( what.str(2) );
        const std::string C( what.str(3) );
        for( size_t site(0); site < 2; ++site )
        { 
          const int
            i = details::type_index()( structure.lattice->sites[ site ? 0: 1 ].type,  A ),
            j = details::type_index()( structure.lattice->sites[ site ].type,  B ),
            k = details::type_index()( structure.lattice->sites[ site ? 0: 1 ].type,  C );
          if( i == -1 or j == -1 or k == -1 ) continue;
          const size_t Bkind
          (
            site == 0 ? size_t(j):( two_species ? size_t(2 + j): size_t(1 + j) )
          );
          const size_t angle_kind
          (
              size_t(site == 0 ? size_t(i):( two_species ? size_t(i): 0 ))
            + size_t(site == 0 ? size_t(k):( two_species ? size_t(k): 0 ))
          );
          LADA_DO_NASSERT( Bkind >= functionals.size(), "Index out-of-range.\n" )
          functionals[Bkind].set_angle( angle_kind, _tuple );
        }
      }

  } // namespace vff 
} // namespace LaDa

#endif // _VFF_FUNCTIONAL_H_
