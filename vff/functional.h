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

#include <tinyxml/tinyxml.h>

#include "lamarck/atom.h"
#include "lamarck/structure.h"
#include "lamarck/lattice.h"
#include "opt/types.h"
#include "opt/opt_function_base.h"

#ifdef _MPI 
  #include "mpi/mpi_object.h"
#endif

//! \brief Reimplements the Valence Force Field %Functional in c++
//! \details Vff, or Valence Force Field functional, is an empirical functional which
//! attempts to model the strain energy of a material from at most three-body
//! interactions. The there body interactions which are considered are
//! bond-stretching, change in bond-angles, and a combination of these two.
//! 
//! Since Vff is, in its present incarnation, at most a three boyd functional,
//! we have built the implementation on a body centered paradigm. In other words,
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

  class Functional;

  //! \brief Represents a single Atom and its first neighbor relationships
  //! \details This class is meant to be used in conjunction with a list of
  //! Ising_CE::Atom, most likely in a Ising_CE::Structure. It contains a
  //! pointer, Atomic_Center::origin, which points a single Ising_CE::Atom. The
  //! first neighbor bonds of this atom are collected as vector of pointers to
  //! Vff::Atomic_Center objects in Atomic_Center::bonds. Since we are concerned
  //! with periodic structures, Atomic_Center::translations and
  //! Atomic_Center::bool record which periodic image of an atom
  //! Atomic_Center::bonds refer to.
  class Atomic_Center
  {
    friend class Functional;
    public:
      class const_iterator;

    protected:
      Ising_CE::Atom *origin; //!< The atom this object is addressing
      //! \brief Other Vff::Atomic_Center objects with which this one is in a bond-relationship
      //! \details Via bonds, a collection of Atomic_Center can be made into a tree,
      //! which can be travelled linearly, or through the first neighbor bonds,
      //! using Atomic_Center::const_iterator.
      //! \sa Atomic_Center::const_iterator, Vff::Functional::construct_centers, 
      //!     Vff::functional::initialize_centers
      std::vector< Atomic_Center* > bonds; 
      //! \brief Allow to relate origin pointer of Atomic_center::bonds to the correct
      //! periodic image with which the bond is made
      std::vector< atat::rVector3d > translations;
      //! \brief A switch to say wether a bond is made directely or with a periodic
      //! image of an Atomic_Center
      std::vector< bool > do_translates;
      //! Ising_CE::Structure on which the Vff  functional is applied.
      Ising_CE :: Structure *structure;
      bool is_site_one; //!< helps determine the kind of atom this is
      bool is_site_one_two_species; //!< helps determine the kind of atom this is
      atat::rVector3d gradient; //!< place holder to compute gradient
      //! atomic index in Ising_CE::Structure::t_Atoms collection of Atomic_Center::structure
      types::t_unsigned index;

    public:
      //! \brief Default Constructor. 
      //! \param _str structure in which \a _e can be found
      //! \param _e atom to which this Atomic_Center relates
      //! \param _i index of _i in _str.atoms collection. Usefull for mpi processing
      Atomic_Center ( Ising_CE::Structure &_str, Ising_CE::Atom &_e, types::t_unsigned _i);
      //! \brief Copy Constructor
      //! \param[in] _c Atomic_Center object to copy
      Atomic_Center   ( const Atomic_Center &_c )
                    : origin(_c.origin), bonds(_c.bonds), translations(_c.translations), 
                      do_translates(_c.do_translates), structure(_c.structure),
                      is_site_one(_c.is_site_one),
                      is_site_one_two_species( _c.is_site_one_two_species), 
                      gradient(0,0,0), index( _c.index) {} 

      //! \brief Returns the kind of atomic center this is
      //! \details In order to use the right coefficients in bond stretching and other
      //! interactions, we have to know the "kind"
      //! functional this is. This will depend on the atomic specie of this
      //! atom and the encoding of Vff::Atomic_Functional array in
      //! Vff::Functional
      types::t_unsigned kind() const;

      //! \brief Adds _bond to  Atomic_Center::bonds if it is in first neighbor relationship
      //! \details This function returns -1 if \a _e is not a bond, and returns the
      //! number of bounds if it is. Note that poeriodic images of _bond are checked and 
      //! Atomic_Center::translations and Atomic_Center::do_translate are set accordingly.
      //! \param _e Atomic_Center object for which first neighbor relationship is checked
      //! \param _cutoff distance below which first neighborness is implied
      types::t_int add_bond( Atomic_Center &_e, const types::t_real _cutoff  );

      //! Returns a Atomic_Center::const_iterator object pointing to the first bond
      const_iterator begin() const;
      //! Returns a Atomic_Center::const_iterator object pointing to the last bond
      const_iterator end() const;
      //! Returns the number of bonds
      types::t_unsigned size() const
        { return bonds.size(); }

      //! Sets the atomic position of the origin
      atat::rVector3d& operator=(atat::rVector3d& _vec)
        { origin->pos = _vec; return _vec; }
      //! Translates the atomic position of the origin by _vec
      //! \param _vec Translation
      void operator+=(const atat::rVector3d& _vec)
        { origin->pos += _vec; }
      //! Translates the atomic position of the origin by -_vec
      //! \param _vec "Negative" Translation
      void operator-=(const atat::rVector3d& _vec)
        { origin->pos -= _vec; }
      //! Returns the atomic position of the origin
      operator atat::rVector3d& ()
        { return origin->pos; }
      //! Returns the atomic position of the origin, constant format
      operator const atat::rVector3d& () const
        { return origin->pos; }
      //! Returns the atom at the origin
      Ising_CE::Atom& get_origin()
        { return *origin; }
      //! Returns the atom at the origin, constant format
      const Ising_CE::Atom& get_origin() const
        { return *origin; }
      //! Returns the gradient place holder
      atat::rVector3d& get_gradient()
        { return gradient; }
      //! Returns the gradient place holder, constant format
      const atat::rVector3d& get_gradient() const
        { return gradient; }
      //! Sets the gradient place-holder to the (0,0,0) vector
      void reset_gradient()
        { gradient[0] = 0; gradient[1] = 0; gradient[2] = 0; }
      //! Returns true if Atomic_Center is a site 1 in this lattice type
      bool site_one() const
        { return is_site_one; }
      //! Returns the index of this Atomic_Center in Atomic_Center::structure
      types::t_unsigned get_index() const
        { return index; }

    protected:
      //! \brief Returns the type of bond between this Atomic_Center and _bond
      //! \param _bond If this is not a bond, function will return result, not error!!
      //! \sa Atomic_Center::add_bond(), Atomic_Functional::add_bond(), Functional::Load()
      types::t_unsigned bond_kind( const Atomic_Center &_bond ) const;
  };

  //! \brief Iterator to travel along the bonds of an Atomic_Center
  //! \details Once a mesh of Atomic_Center objects is constructed, with the bonds
  //! forming the links in the mesh, an iterator is needed which will allow us to
  //! travel the mesh in any direction. Atomic_Center::const_iterator is built
  //! to iterates through the bonds of an Atomic_Center object. From there, one
  //! can easily travel throughout the mesh via first neighbor relationships.
  class Atomic_Center :: const_iterator
  {
    protected:
      //! current origin of the bonds
      const Atomic_Center *parent;
      //! current bond being iterated
      std::vector< Atomic_Center* > :: const_iterator i_bond;
      //! tranlation, if current bond is an atomic image
      std::vector< atat::rVector3d > :: const_iterator i_translation;
      //! swtich for doing periodic imeage translation or not
      std::vector< bool > :: const_iterator i_do_translate;
    
    public:
      //! Constructor.
      const_iterator() {};
      //! \brief Constrcutor and Initialized
      //! \param _parent origin of the bond
      //! \param _is_begin
      //!                   - on true, initializes to first bond
      //!                   - on false, initializes to last bond
      //!                   .
      const_iterator   ( const Atomic_Center *_parent, bool _is_begin=true) 
                     : parent(_parent)
      {
        if ( _is_begin )
        {
          i_bond = parent->bonds.begin();
          i_translation = parent->translations.begin();
          i_do_translate = parent->do_translates.begin();
          return;
        }
        i_bond = parent->bonds.end();
      }
      //! Copy Constrcutor 
      const_iterator   ( const const_iterator &_c ) 
                     : parent(_c.parent), i_bond(_c.i_bond), 
                       i_translation(_c.i_translation),
                       i_do_translate(_c.i_do_translate) {};
      //! Iterates through bonds
      void operator ++()
        { ++i_bond; ++i_translation; ++i_do_translate; }
      //! Iterates through bonds
      void operator --()
        { --i_bond; --i_translation; --i_do_translate; }
      //! Iterates through bonds
      void operator -=( types::t_int  n)
        { i_bond -= n; i_translation -= n; i_do_translate -= n; }
      //! Iterates through bonds
      void operator +=( types::t_int  n)
        { i_bond += n; i_translation += n; i_do_translate += n; }
      //! \brief Returns true if iterators are at same bond
      //! \param _i iterator agains which to check
      bool operator ==( const const_iterator &_i ) const
        { return _i.i_bond == i_bond; }
      //! \brief Returns true if iterators are \em not at same bond
      //! \param _i iterator agains which to check
      bool operator !=( const const_iterator &_i ) const
        { return _i.i_bond != i_bond; }
      //! \brief Returns the number of bonds separating this object  and _i in
      //!  bond array 
      //! \param _i iterator agains which to check
      types::t_int operator -( const const_iterator &_i )
        { return _i.i_bond - i_bond; }
      //! Dereferences the iterator: returns Atomic_Center to which bond points
      Atomic_Center& operator *()
        { return *(*i_bond); }
      //! Dereferences the iterator: returns Atomic_Center to which bond points
      const Atomic_Center& operator *() const
        { return *(*i_bond); }
      //! Dereferences the iterator: returns Atomic_Center to which bond points
      Atomic_Center* operator ->()
        { return *i_bond; }
      //! Dereferences the iterator: returns Atomic_Center to which bond points
      const Atomic_Center* operator ->() const
        { return *i_bond; }
      //! Returns bond length squared
      types::t_real norm2() const
      {
        if ( not *i_do_translate )
          return atat::norm2 ( parent->origin->pos - (*i_bond)->origin->pos );
        return atat::norm2( parent->origin->pos - (*i_bond)->origin->pos -
                            parent->structure->cell * (*i_translation) );
      }
      //! \brief Returns bond vector
      atat::rVector3d& vector( atat::rVector3d &_hold )
      {
        _hold = (*i_bond)->origin->pos - parent->origin->pos ;
        if ( *i_do_translate )
          _hold += parent->structure->cell * (*i_translation);
        return _hold;
      }
      //! \brief Returns scalar product between bond vector and _b
      types::t_real scalar_product( const const_iterator &_b ) const
      {
        atat::rVector3d a, b;
        if ( *i_do_translate )
          a =   (*i_bond)->origin->pos - parent->origin->pos 
              + parent->structure->cell * (*i_translation);
        else
          a = (*i_bond)->origin->pos - parent->origin->pos;
        if ( *_b.i_do_translate )
          b =   (*_b.i_bond)->origin->pos - _b.parent->origin->pos 
              + _b.parent->structure->cell * (*_b.i_translation);
        else
          b = (*_b.i_bond)->origin->pos - _b.parent->origin->pos;
        return a * b;
      }
      //! \brief Returns the kind of bond this is
      //! \see  Atomic_Center::bond_kind(), Atomic_Center::add_bond(),
      //!       Atomic_Functional::add_bond(), Functional::Load() 
      types::t_unsigned kind() const
        { return parent->bond_kind( *(*i_bond) ); }
      //! \brief Returns the atom at the origin of range of bonds this iterator travels
      //! \see Atomic_Center::const_iterator::parent 
      Ising_CE::Atom& get_origin()
        { return ((*i_bond)->get_origin()); }
      //! \brief Translates a vector _v by periodic image of enpoint of bond
      //! \param _v vector to translate
      //! \param _cell unit-cell defining periodic image (can be different from
      //! Atomic_Center::structure by amount of minimized strain)
      void translate( atat::rVector3d &_v, const atat::rMatrix3d &_cell )
        { if( *i_do_translate ) _v += _cell * ( *i_translation ); }
      //! \brief Translates a vector _v by periodic image of enpoint of bond
      //! \param _v vector to translate
      void translate( atat::rVector3d &_v )
        { if( *i_do_translate ) _v += parent->structure->cell * ( *i_translation ); }
  }; // end of const_iterator definition

  //! \brief An atom centered functional
  //! \details Atomic_Functional will compute the energy and strain resulting from
  //! bond-stretching, angle warping, and bond-angle-warping, eg from all
  //! first-neighbor two and three body interactions, on an Vff:;Atomic_Center atom.
  class Atomic_Functional 
  {
#ifdef _MPI
    /// \cond
    friend bool mpi::BroadCast::serialize<Vff::Atomic_Functional> ( Vff::Atomic_Functional& );
    /// \endcond
#endif
     const static types::t_real twos3;  //!<\f$2*\sqrt(3)\f$
     const static types::t_real one16;  //!<\f$\frac{1}{16}\f$
     const static types::t_real s3o160; //!<\f$\frac{\sqrt(3)}{8}\f$
     const static types::t_real one640; //!<\f$\frac{1}{640}\f$
     const static types::t_real three8; //!<\f$\frac{3}{8}\f$
     const static types::t_real s33o8;  //!<\f$\frac{3}{8}\sqrt(3)\f$
     const static types::t_real s33o16; //!<\f$\frac{3}{16}\sqrt(3)\f$
     const static types::t_real thre16; //!<\f$\frac{3}{16}\f$
     const static types::t_real thre32; //!<\f$\frac{3}{32}\f$
     const static types::t_real s33128; //!<\f$\frac{3}{128}\sqrt(3)\f$
     const static types::t_real s33256; //!<\f$\frac{3}{256}\sqrt(3)\f$
     const static types::t_real no1280; //!< some number, but which?
     const static types::t_real no2560; //!< some number, but which?

    protected:
      std::string str;                  //!< atomic type as a string
      Ising_CE :: Structure *structure; //!< structure to which the Atomic_Center belongs
      types::t_unsigned site;           //!< site number of Atomic_Center in Ising_CE::lattice
      types::t_unsigned type;           //!< atomic type of Atomic_Center
      std::vector< types::t_real > lengths;  //!< equilibrium bond lengths
      std::vector< types::t_real > alphas;   //!< bond stretching parameters
      std::vector< types::t_real > betas;    //!< angle deformation parameters
      std::vector< types::t_real > gammas;   //!< bond-angle parameters
      std::vector< types::t_real > sigmas;   //!< equilibrium tetrahedral symmetry
      
    public:
#ifdef _MPI
      /** \ingroup MPI
      *   \brief Constructor,  should only be used by serialize
      *   \details Can't remember why anymore...
      */ 
      Atomic_Functional   ( Ising_CE::Structure &_struct ) 
                        : structure(&_struct) {}  
#endif
      //! \brief Constructor and Initializer
      //! \param _str atomic type as a string
      //! \param _struct structure to which Atomic_Center belongs
      //! \param _site site number of Atomic_Center
      //! \param _type type number of Atomic_Center
      Atomic_Functional   ( std::string _str, Ising_CE::Structure &_struct, 
                            types::t_unsigned _site, 
                            types::t_unsigned _type )
                        : structure(&_struct), site(_site), type(_type) {}
      //! Copy Constructor
      Atomic_Functional   ( const Atomic_Functional &_a )
                        : str(_a.str), structure(_a.structure), site(_a.site), type(_a.type),
                          lengths(_a.lengths), alphas(_a.alphas),
                          betas(_a.betas), gammas(_a.gammas), sigmas(_a.sigmas) {}
      
      //! \brief Adds a bond type to bond list
      //! \param _typeB type of atom at end-point of bond
      //! \param _l equilibrium length
      //! \param _i bond stretching parameters
      void add_bond( const types::t_unsigned _typeB, const types::t_real _l,
                     const std::vector<types::t_real> &_i )
      {
        if ( lengths.size() < _typeB + 1)
          lengths.resize( _typeB+1, types::t_real(0) );
        if ( alphas.size() < (_typeB+1) * 5  )
          alphas.resize( (_typeB+1)*5, types::t_real(0) );
        lengths[_typeB] = _l;
        std::copy( _i.begin(), _i.end(), alphas.begin() + _typeB * 5 );
      }
      //! \brief Adds angle deformation and bond-angle parameters
      //! \param _typeA type of atom at end-point of one bond
      //! \param _typeC type of atom at end-point of other bond
      //! \param _gamma bond-angle deformation parameter
      //! \param _sigma equilibrium angle
      //! \param _i angle deformation parameters
      void add_angle( const types::t_unsigned _typeA,
                      const types::t_unsigned _typeC,
                      const types::t_real _gamma, const types::t_real _sigma, 
                      const std::vector<types::t_real> &_i )
      {
        types::t_unsigned offset = _typeA+_typeC;
        if ( gammas.size() < offset + 1  )
          gammas.resize( offset + 1, types::t_real(0) );
        if ( sigmas.size() < offset + 1  )
          sigmas.resize( offset + 1, types::t_real(0) );
        if ( betas.size() < (offset + 1) * 5 )
          betas.resize( (offset+1)*5, types::t_real(0) );
        gammas[offset] = _gamma;
        sigmas[offset] = _sigma;
        std::copy( _i.begin(), _i.end(), betas.begin() + offset * 5 );
      }
      
      //! \brief Evaluate strain energy for Atomic_Center _center
      //! \details returns strain energy 
      //! \param _center center for which to evaluate energy
      //! \sa function::Base, function::Base::evaluate()
      types::t_real evaluate( const Atomic_Center &_center ) const;
      //! \brief Evaluate strain energy and gradients for Atomic_Center _center
      //! \details returns strain energy, and computes stress
      //! \param _center center for which to evaluate energy
      //! \param _strain to which Atomic_Functional::structure is submitted
      //! \param _stress on _center resulting from _strain
      //! \param _K0 displacement resulting from _strain and w.r.t original unit-cell
      //! \sa function::Base, function::Base::evaluate_with_gradient()
      types::t_real evaluate_with_gradient( Atomic_Center &_center,
                                            const atat::rMatrix3d &_strain,
                                            atat::rMatrix3d &_stress,
                                            const atat::rMatrix3d &_K0 ) const;
      //! \brief computes the trace of the microscopic strain on an atomic center
      //!  structure0 and the atomic centers are expected to be related 
      //! \details To be used for pescan
      types::t_real MicroStrain( const Atomic_Center &_center, 
                                 const Ising_CE::Structure &_str0 ) const;

      //! prints out all parameters
      void print_out( std::ostream &stream ) const;
  }; 

  //! Wrapper, or interface class to Valence Force Field Functional
  //
  //! Takes care of loading parameters from input, building tree for Functional::structure,
  //! and eventually is able to compute its the energy and strain.
  //! One possible use is the following:
  //! \code
  //! Ising_CE::Structure structure;
  //! Ising_CE::Lattice lattice;
  //! Ising_CE :: Structure :: lattice = &lattice; // don't forget this hack
  //! // load lattice, structure, etc from input
  //! // Then load Vff input itself
  //! TiXmlElement *vff_xml = handle.FirstChild( "Job" ).Element();
  //! Vff::Functional vff(structure);
  //! if ( not vff.Load(*vff_xml) ) return false;
  //! 
  //! vff.initialize_centers(); // Construct mesh of first neighbor relationships
  //!
  //! 
  //! minimizer::GnuSL<Vff::Functional> minimizer( vff ); // Declare and Load a GSL minimizer
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
#ifdef _MPI
    /// \cond
    friend bool mpi::BroadCast::serialize<Vff::Functional> ( Vff::Functional& );
    /// \endcond
#endif
    public:
      typedef types::t_real t_Type;            //!< see Functional::Base
      typedef std::vector<t_Type> t_Container; //!< see Functional::Base
      typedef t_Container :: iterator iterator; //!< see Functional::Base
      typedef t_Container :: const_iterator const_iterator; //!< see Functional::Base

    protected:
      //! Ising_CE::Structure for which to compute energy and stress
      Ising_CE :: Structure &structure;
      Ising_CE :: Structure structure0; //!< original structure,  needed for gradients
      types::t_real bond_cutoff; //!< length below which first-neighbor relationship is defined
      //! list of all Atomic_Center created from Functional::structure
      std::vector< Atomic_Center > centers;  
      //! list of all possbile Atomic_Functionals for Functional::structure.lattice
      std::vector< Atomic_Functional > functionals;
      atat::rMatrix3d stress; //!< stores stress in Functional::structure after computation
      atat::rMatrix3d strain; //!< stores stress in Functional::structure after computation
      atat::rVector3d center_of_mass; //!< center of mass of atoms in Functional::structure
      
    public:
      //! \brief Constructor and Initializer
      //! \param _str structure for which to compute energy and stress
      Functional   ( Ising_CE :: Structure &_str )
                 : function::Base<>( 7 + _str.atoms.size() ),
                   structure(_str), structure0(_str), center_of_mass(0,0,0) {};
      //! \brief Copy Constructor
      Functional   ( const Vff::Functional &_c )
                 : function::Base<>( _c ), structure( _c.structure ),
                   structure0( _c.structure0 ), bond_cutoff( _c.bond_cutoff ),
                   centers( _c.centers ), functionals( _c.functionals ),
                   center_of_mass(_c.center_of_mass) {}
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
      //! \sa function::Base, function::Base::evaluate_one_gradient, minimizer::VA
      t_Type evaluate_one_gradient( types::t_unsigned _pos) {return 0;}; 
      //! \brief initializes stuff before minimization
      //! \details Defines the packing and unpacking process, such that only unfrozen
      //! degrees of liberty are known to the minimizer
      //! \sa function::Base, minimizer::Base
      bool init();
      //! \brief Prints atom.config type input to escan
      //! \param _f optional filename to which to direct the output
      void print_escan_input( const std::string &_f = "atom.config") const;
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
      //! \brief packs variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas minizers now
      //! about function::Base, this function does the interface between the two
      void pack_variables(const atat::rMatrix3d& _strain);
      //! \brief packs variables from minimizer
      //! \details Functional knows about Functional::Structure, whereas
      //! minizers now about function::Base, this function does the interface
      //! between the two
      template< typename t_grad_iterator>
      void pack_gradients(const atat::rMatrix3d& _stress, t_grad_iterator const &_grad) const;
  };

  template< typename t_grad_iterator>
  void Functional :: pack_gradients(const atat::rMatrix3d& _stress, 
                                    t_grad_iterator const &_grad) const
  {
    t_grad_iterator i_grad(_grad);

    // first, external stuff
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_XX) )
      *i_grad = _stress(0,0), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_YY) ) 
      *i_grad = _stress(1,1), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_ZZ) ) 
      *i_grad = _stress(2,2), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_XY) ) 
      *i_grad = 0.5 * (_stress(0,1) + _stress(1,0)), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_XZ) ) 
      *i_grad = 0.5 * (_stress(0,2) + _stress(2,0)), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_YZ) ) 
      *i_grad = 0.5 * (_stress(1,2) + _stress(2,1)), ++i_grad;

    // then atomic position stuff
    std::vector<Atomic_Center> :: const_iterator i_center = centers.begin();
    std::vector<Atomic_Center> :: const_iterator i_end = centers.end();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom0 = structure0.atoms.begin();
    i_center = centers.begin();
    for (; i_center != i_end; ++i_center, ++i_atom0)
    {
      const atat::rVector3d& gradient = i_center->get_gradient();
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_X) ) 
        *i_grad = gradient[0], ++i_grad;
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_Y) ) 
        *i_grad = gradient[1], ++i_grad;
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_Z) ) 
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
    std::vector<Atomic_Center> :: iterator i_center = centers.begin();
    std::vector<Atomic_Center> :: iterator i_end = centers.end();
    stress.zero();
    for (; i_center != i_end; ++i_center)
      energy += functionals[i_center->kind()].evaluate_with_gradient( *i_center, strain, stress, K0 );

    // now repacks into function::Base format
    pack_gradients(stress, _i_grad);

    return energy;
  }

} // namespace vff 

namespace mpi {
  /** \ingroup MPI
  * \brief Serializes Vff::Atom_Functional. Not sure this works anymore.
  */
  template<>
  bool BroadCast::serialize<Vff::Atomic_Functional> ( Vff::Atomic_Functional& );
  
  /** \ingroup MPI
  * \brief Serializes Vff::Atom_Functional, meaning all member variables
  * without exception.
  */
  template<>
  bool BroadCast::serialize<Vff::Atomic_Functional> ( Vff::Atomic_Functional& );
}
#endif // _VFF_FUNCTIONAL_H_
