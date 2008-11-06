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

//! \brief Reimplements the Valence Force Field %Functional in c++
//! \details Vff, or Valence Force Field functional, is an empirical functional which
//! attempts to model the strain energy of a material from at most three-body
//! interactions. The there body interactions which are considered are
//! bond-stretching, change in bond-angles, and a combination of these two.
//! 
//! The implementation relies on a body-centered paradigm. In other words,
//! four classes have been created:
//!   - Vff::Functional is wrapper class and interface to Vff
//!   - Vff::Atomic_Center represent a single atom and lists its first neighbor relationships
//!   - Vff::Atomic_Center::const_iterator allows coders to travel
//!   through a collection of Vff::Atomic_Centera along first neighbor
//!   relationships.
//!   - Vff::Atomic_Functional computes the strain energy of one single atom, eg
//!   all three body terms in which a particular atom takes part.
//!   .
namespace Vff
{
  //! \brief Wrapper, or interface class to Valence Force Field Functional
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
  class Functional : public function :: Base<types::t_real, std::vector<types::t_real> >
  {
    //! Type of the path.
    typedef boost::filesystem::path t_Path;
    //! The type of the atom  
    typedef Crystal::Structure::t_Atom  t_Atom;
    //! The type of the atom container
    typedef Crystal::Structure::t_Atoms t_Atoms;
    public:
      typedef types::t_real t_Type;            //!< see Functional::Base
      typedef std::vector<t_Type> t_Container; //!< see Functional::Base
      typedef t_Container :: iterator iterator; //!< see Functional::Base
      typedef t_Container :: const_iterator const_iterator; //!< see Functional::Base

    protected:
      //! Type of the atomic centers
      typedef Atomic_Center t_Center;  
      //! Type of the container holding the atomic centers
      typedef t_Center :: t_Centers t_Centers;  
      //! Type of the container holding the atomic functionals
      typedef std::vector< Atomic_Functional > t_AtomicFunctionals;  


    protected:
      //! Crystal::Structure for which to compute energy and stress
      Crystal :: Structure &structure;
      Crystal :: Structure structure0; //!< original structure,  needed for gradients
      //! length below which first-neighbor relationship is defined
      types::t_real bond_cutoff; 
      //! \brief list of all Atomic_Center created from Functional::structure
      //! \details Space for the atomic centers are reserved in the
      //!          constructor. It is expected that the iterators will be valid
      //!          throughout the life of the functional, starting with a call
      //!          to the construction of the tree. If the order of the atoms
      //!          in Functional::structure is changed, then it is imperative
      //!          that the tree be reconstructed from scratch.
      t_Centers centers;  
      //! list of all possbile Atomic_Functionals for Functional::structure.lattice
      t_AtomicFunctionals functionals;
      //! stores stress in Functional::structure after computation
      atat::rMatrix3d stress;
      //! stores stress in Functional::structure after computation
      atat::rMatrix3d strain; 
      //! Index of the first atoms with fixed x, y, z;
      atat::iVector3d fixed_index; 
      
    public:
      //! \brief Constructor and Initializer
      //! \param _str structure for which to compute energy and stress
      Functional   ( Crystal :: Structure &_str )
                 : structure(_str), structure0(_str),
                   bond_cutoff(0), fixed_index(-1,-1,-1)
      {
        functionals.reserve( _str.atoms.size()); 
        stress.zero(); strain.zero();
      };
      //! \brief Copy Constructor
      Functional   ( const Vff::Functional &_c )
                 : function::Base<>( _c ), structure( _c.structure ),
                   structure0( _c.structure0 ), bond_cutoff( _c.bond_cutoff ),
                   centers( _c.centers ), functionals( _c.functionals ),
                   fixed_index( _c.fixed_index ) {}
      //! \brief Destructor
      ~Functional() {}

      //! \brief Loads input to functional from  xml 
      //! \param _element should point to an xml node which is the functional data
      //! or contains a node to the funtional data as a direct child
      bool Load( const TiXmlElement &_element );

      //! \brief unpacks function::Base::variables, then calls energy
      //! \sa function::Base, function::Base::evaluate
      types::t_real evaluate(); 
      //! \brief computes energy and stress, expects everything to be set
      types::t_real energy() const;
      //! \brief Evaluates gradients only
      //! \sa function::Base, function::Base::evaluate_gradient
      template< typename t_grad_iterator>
        void evaluate_gradient( t_grad_iterator const &_i_grad )
          { evaluate_with_gradient( _i_grad ); }
      //! \brief Evaluates gradients only
      //! \sa function::Base, function::Base::evaluate_gradient
      void evaluate_gradient( t_Type * const _i_grad )
        { evaluate_with_gradient<t_Type*>( _i_grad ); }  
      //! \brief Evaluates gradients and energy
      //! \sa function::Base, function::Base::evaluate_with_gradient
      template< typename t_grad_iterator>
        t_Type evaluate_with_gradient( t_grad_iterator const &_i_grad );
      //! \brief Evaluates gradients and energy
      //! \sa function::Base, function::Base::evaluate_with_gradient
      t_Type evaluate_with_gradient( t_Type * const _i_grad )
        { return evaluate_with_gradient<t_Type*>( _i_grad ); }  
      //! \brief Evaluates gradient in one direction only
      //! \todo Vff::Functional::implement evaluate_one_gradient
      //! \sa function::Base, function::Base::evaluate_one_gradient, Minimizer::VA
      t_Type evaluate_one_gradient( types::t_unsigned _pos) {return 0;}; 
      //! \brief initializes stuff before minimization
      //! \details Defines the packing and unpacking process, such that only unfrozen
      //! degrees of liberty are known to the minimizer
      //! \sa function::Base, Minimizer::Base
      bool init();
      //! \brief Prints atom.config type input to escan
      //! \param _f optional filename to which to direct the output
      void print_escan_input( const t_Path &_f = "atom.config") const;
      //! \brief Constructs the mesh of Atomic_Center
      //! \details Newer version than Functional::initialize_centers, it works even if
      //! structure is slightly distorted.
      bool construct_centers();
      //! \deprecated Constructs the mesh of Atomic_Center
      bool initialize_centers();
      //! Prints out all parameters
      void print_out( std::ostream &stream ) const;
      //! \brief Returns a reference to the computed stress
      //! \sa Functional::stress
      const atat::rMatrix3d& get_stress() const
        { return stress; }

    protected:
      //! \brief unpacks variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas minizers now
      //! about function::Base, this function does the interface between the two
      void unpack_variables(atat::rMatrix3d& strain);
      //! Unpacks position variables only.
      void unpack_positions(atat::rMatrix3d& strain, const_iterator& _i_x);
      //! \brief packs variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas minizers now
      //! about function::Base, this function does the interface between the two
      void pack_variables(const atat::rMatrix3d& _strain);
      //! Packs position variables only.
      void pack_positions(iterator& _i_x);
      //! Counts positional degrees of freedom.
      types::t_unsigned posdofs();
      //! \brief packs variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas
      //! minizers now about function::Base, this function does the interface
      //! between the two
      template< typename t_grad_iterator>
      void pack_gradients(const atat::rMatrix3d& _stress,
                          t_grad_iterator const &_grad) const;
      //! \brief Loads Functional for one or two site lattices.
      //! \details If \a _node is not the correct node, the results are undefined.
      bool load_( const TiXmlElement &_node );
      

#ifdef _LADADEBUG
      //! Checks that the list of centers are valid. somewhat.
      void check_tree() const;
#endif
  };

  template< typename t_grad_iterator>
  void Functional :: pack_gradients(const atat::rMatrix3d& _stress, 
                                    t_grad_iterator const &_grad) const
  {
    t_grad_iterator i_grad(_grad);

    // first, external stuff
    if ( not (structure.freeze & Crystal::Structure::FREEZE_XX) )
      *i_grad = _stress(0,0), ++i_grad;
    if ( not (structure.freeze & Crystal::Structure::FREEZE_YY) ) 
      *i_grad = _stress(1,1), ++i_grad;
    if ( not (structure.freeze & Crystal::Structure::FREEZE_ZZ) ) 
      *i_grad = _stress(2,2), ++i_grad;
    if ( not (structure.freeze & Crystal::Structure::FREEZE_XY) ) 
      *i_grad = 0.5 * (_stress(0,1) + _stress(1,0)), ++i_grad;
    if ( not (structure.freeze & Crystal::Structure::FREEZE_XZ) ) 
      *i_grad = 0.5 * (_stress(0,2) + _stress(2,0)), ++i_grad;
    if ( not (structure.freeze & Crystal::Structure::FREEZE_YZ) ) 
      *i_grad = 0.5 * (_stress(1,2) + _stress(2,1)), ++i_grad;

    // then atomic position stuff
    t_Centers :: const_iterator i_center = centers.begin();
    t_Centers :: const_iterator i_end = centers.end();
    t_Atoms :: const_iterator i_atom0 = structure0.atoms.begin();
    i_center = centers.begin();
    for (; i_center != i_end; ++i_center, ++i_atom0)
    {
      const atat::rVector3d& gradient = i_center->get_gradient();
      if ( not (i_atom0->freeze & t_Atom::FREEZE_X) ) 
        *i_grad = gradient[0], ++i_grad;
      if ( not (i_atom0->freeze & t_Atom::FREEZE_Y) ) 
        *i_grad = gradient[1], ++i_grad;
      if ( not (i_atom0->freeze & t_Atom::FREEZE_Z) ) 
        *i_grad = gradient[2], ++i_grad;
    }
  }

  template< typename t_grad_iterator>
  types::t_real Functional :: evaluate_with_gradient( t_grad_iterator const &_i_grad )
  {
    t_Type energy = 0;
    std::for_each( centers.begin(), centers.end(), std::mem_fun_ref(&Atomic_Center::reset_gradient) );

    // unpacks variables into vff atomic_center and strain format
    unpack_variables(strain);

    // computes K0
    atat::rMatrix3d K0 = (!(~strain));

    // computes energy and gradient
    t_Centers :: iterator i_center = centers.begin();
    t_Centers :: iterator i_end = centers.end();
    stress.zero();
    for (; i_center != i_end; ++i_center)
      energy += functionals[i_center->kind()].
                     evaluate_with_gradient( *i_center, strain, stress, K0 );

    // now repacks into function::Base format
    pack_gradients(stress, _i_grad);

    return energy;
  }

  // same as energy, but unpacks values from
  // opt::Function_Base::variables
  inline types::t_real Functional :: evaluate()
  {
    unpack_variables(strain);
    return energy();
  }
} // namespace vff 

#ifdef _DOFORTRAN
#include<opt/opt_frprmn.h>
  //! Creates an instance of a typical Minimizer::Frpr "C" function for calling Vff::Functional
  extern "C" inline double vff_frprfun(double* _x, double* _y)
    { return Minimizer::typical_frprfun<Vff::Functional>( _x, _y); }
  //! \brief returns a pointer to the correct extern "C" evaluation function
  //!        for Minimizer::Frpr.
  //! \details This routine allows for a standard for Vff::VA to intialize
  //!          Minimizer::Frpr.
  template<> inline t_FrprFunction choose_frpr_function<Vff::Functional>() { return vff_frprfun; }

#elif defined(_DONAG)
#include <nag.h>
#include <nage04.h>
  //! Creates an instance of a typical NAG "C" function for calling Vff::Functional
  extern "C" inline void vff_nagfun(int _n, double* _x, double* _r, double* _g, Nag_Comm* _p)
    { Minimizer::typical_nagfun<Vff::Functional>( _n, _x, _r, _g, _p ); }
  //! \brief returns a pointer to the correct extern "C" evaluation function
  //!        for Minimizer::Nag.
  //! \details This routine allows for a standard for Vff::VA to intialize
  //!          Minimizer::Nag.
  template<>
  inline Minimizer::t_NagFunction choose_nag_function<Vff::Functional>()
    { return vff_nagfun; }

#endif

#endif // _VFF_FUNCTIONAL_H_
