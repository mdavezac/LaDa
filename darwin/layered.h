//
//  Version: $Id$
//
#ifndef _LAYERED_H_
#define _LAYERED_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <ostream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <lamarck/structure.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "evaluator.h"
#include "individual.h"
#include "bitstring.h"
#include "gaoperators.h"
#include "vff.h"


//! \brief fills in \a atoms member of an Ising_CE::Structure instance from the
//!        cell-shape and the lattice.
//! \todo make sure that the guessed range is large enough, so that any
//!       structure is completely filled with atoms.
void FillStructure( Ising_CE::Structure &_str );

/** \ingroup Genetic
 * @{ */
//! \brief allows the optimization of a (specified) structure's decoration,
//!        where the structure \e layered.
//! \details The point is to optimize for epitaxially grown structures. As such
//!          a structure is defined simply by:
//!             - The lattice type 
//!             - An epitaxial growth direction.
//!             . The number of independent layers
//!             .
//!          From this information a complete instance of anIsing_CE::Structure
//!          is constructed, The atoms are ordered according to depth in the
//!          epitaxial growth direction. The structure can then be mapped onto
//!          a bitstring in which bits correspond to atomic-site with
//!          (decreasing?) increasing depth along the epitaxial growth
//!          direction.
//!
//!          Furthermore, in the case of quaternaries, layers are grown
//!          coherently. Two one cation correspond one anion. To the other
//!          cation corresponds the other anion. The bitstring of
//!          Layered::Object<2> contains only half as many atoms as the
//!          Ising_CE::Structure which it characterizes. Note that all sites -
//!          other than the first - are set to have  their occupations frozen.
//!          This allows us to work correctly with the physical operators
//!          defined in gaoperators.h.
//!          
//!          In layered are defined the fourier transform, concentration,
//!          object, and evaluator stub classes capable of dealing with layered
//!          structures for any epitaxial growth direction, any lattice, with
//!          any number of lattice sites. (a layer is always composed of a
//!          complete complement of lattice-sites.) 
namespace Layered
{
  //! Defines %Fourier transforms for a layered structure of a lattice with \a _D sites
  template<types::t_unsigned _D>
  struct Fourier
  {
    //! \brief From real to k space
    //! \details The first range is most likely some instance of
    //!          Ising_CE::Structure::t_Atoms. The second range is similarly an
    //!          instance of Ising_CE::Structure::t_kAtoms. The third range
    //! \param[in, out] _rfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type)
    //! \param[in, out] _rend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type< std::complex >)
    template<class T_R_IT, class T_K_IT>
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend );
    //! \brief From k space to real space. 
    //! \details The first range is most likely some instance of
    //!          Ising_CE::Structure::t_Atoms. The second range is similarly an
    //!          instance of Ising_CE::Structure::t_kAtoms. The third range
    //!          should be iterators to std::complex.
    //! \pre The range [ \a _rout, \a _rout += \a _rfirst - \a _rend  ) should
    //!      be valid.
    //! \param[in] _rfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type)
    //! \param[in] _rend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Ising_CE::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Ising_CE::Atom_Type< std::complex >)
    //! \param[out] _rout iterator to the first complex real-space
    //!              occupation ( of std::complex type )
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend,
             T_O_IT _rout ); // sets rvector values from kspace values
  };


  //! \brief %Layered %Object Type
  //! \details Redefines BitString::Object with the sole purpose of implementing
  //!          overloaded operator<<() on and from Ising_CE::Structures.
  //!          Note that we could run into awkward redefinitions by using typedefs
  //!          rather than a complete redeclaration.
  //! \see Layered::operator<<( Ising_CE::Structure&, const Layered::Object& ), 
  //!      Layered::operator<<( Layered::Object&, const Ising_CE::Structure& ). 
  template<class T_CONTAINER = std::vector<types::t_real> >
  class Object: public BitString :: Object<T_CONTAINER>
  {
    public:
      //! \brief Type of container for which to do bitstring %GA
      //! \details See also function::Base::t_Container
      typedef T_CONTAINER t_Container;
      //! \brief See function::Base::t_Type
      typedef typename t_Container :: value_type t_Type;
      //! \brief See function::Base::iterator
      typedef typename t_Container :: iterator iterator;
      //! \brief See function::Base::const_iterator
      typedef typename t_Container :: const_iterator const_iterator;
    protected:
      //! Type of the base class 
      typedef BitString :: Object<T_CONTAINER> t_Base;

    public:
      //! Constructor
      Object() : t_Base() {}
      //! Copy Constructor
      Object(const Object &_c) : t_Base(_c) {};
      //! Constructor and Initializer
      Object(const t_Container &_c) : t_Base(_c) {};
      //! Destructor
      ~Object() {};
  };

  //! Dumps the decoration of \a _str into the object \a _o
  template<class T_CONTAINER>
    void operator<<(Ising_CE::Structure &_str, const Object<T_CONTAINER> &_o);
  //! Dumps the decoration of \a _str into the object \a _o
  template<class T_CONTAINER>
    void operator<<(Object<T_CONTAINER> &_o, const Ising_CE::Structure &_str);

  //! \brief Declares concentration related behaviors for layered structures.
  //! \details It can compute and set the concentrations of Object and
  //!          Ising_CE::Structure instances for layered structures. It works
  //!          correctly with structure for which the occupation of some sites
  //!          are frozen.
  //! \xmlinput Specifying an x0 attribute in the \<%GA\> tag implies that the
  //!           concetration will be set to that value throughout.
  template <types::t_unsigned _D>
  class Concentration
  {
    public:
      //! Number of sites in this lattice
      const types::t_unsigned _d;

    protected:
      //! \brief Value of the fixed concentration.
      //! \details Not used if the concentration is not fized
      types::t_real x0;
      //! The number of sites in the current cell-shape
      types::t_unsigned N;
      //! \brief The number of sites \e with \e fixed \e occupations in the
      //!        current cell-shape
      types::t_int Nfreeze;
      //! True if the concentration is fixed
      bool single_c;

    public:
      //! \brief Stores the concentation of the last computed Object of
      //!        Ising_CE::Structure instance..
      types::t_real x;

    public:
      //! Constructor
      Concentration() : _d(_D), x0(0.0), N(0), Nfreeze(0),
                        single_c(false), x(0) {}
      //! Copy Constructor
      Concentration   ( const Concentration &_c)
                    : _d(_D), x0(_c.x0), N(_c.N), Nfreeze(_c.Nfreeze),
                      single_c(_c.single_c), x(_c.x) {}
      //! Destructor
      ~Concentration() {}


      //! \brief Normalizes the site occupations as given by the \b k-vectors. 
      //! \details A "complex" real-space occupation is computed from the
      //!          k-space intensities. Normalized site-occupations are set to
      //!          +/- 1 depending on which half-plane the complex value are.
      //!          If the concentration is fixed, those sites for which the
      //!          real value of the complex occupation  are closest to zero
      //!          are flipped first.
      //! \note The concentration is computed from the occupation of the first
      //!       lattice-site only. The occupation of other lattice sites depend
      //!       upon the occupation of the corresponding first lattice-site
      //! \see GA::Krossover, GA::KRandom, GA::KMutation.
      void operator()( Ising_CE::Structure &_str );
      //! \brief Sets \a _obj to concentration Concentration::x0 if the concentration
      //!        is fixed.
      template<class T_CONT> void operator()( BitString::Object<T_CONT> &_obj );
      //! \brief Computes the concentration of \a _str and stores the result in
      //!        Concentration::x
      void get( const Ising_CE::Structure &_str);
      //! \brief Computes the concentration of \a _obj and stores the result in
      //!        Concentration::x
      template<class T_CONT> void get( const BitString::Object<T_CONT> &_obj );
      //! \brief Computes the number of sites in \a _str for which the
      //!        occupation is fixed
      void setfrozen ( const Ising_CE::Structure &_str );

      //! \brief Returns a string describing the statuc of this instance
      std::string print() const;

      //! \brief Chechs whether an attribute is x0. If it is, the concentration
      //!        is fixed.
      //! \details This is used mainly to load attributes from the \<GA\> Tag
      //!          by GA::Darwin.
      void LoadAttribute ( const TiXmlAttribute &_att );

      //! Returns true if the concentration is fixed 
      bool is_single_c() const { return single_c; }


    protected:
      //! \brief Actually does the job of normalizing the site occupations for
      //!        Concentration::operator()(Ising_CE::Structure&)
      void normalize( Ising_CE::Structure &_str, 
                      types::t_real _tochange);

  };

  //! \brief Partially overrides GA::Evaluator default to implement behaviors
  //!        relevant to layered, epitaxial structure decoration search. 
  //! \details Mostly allows for physical %GA operators (GA::Krossover,
  //!          GA::xTaboo, GA::Crossover, etc...) and links them with the
  //!          bitstring objects defined above, as well as the concentration
  //!          functor SingleSite::Concentration. Saving and Restarting of
  //!          individuals is partially implementated as relates to
  //!          BitString::Object and Ising_CE::Structure.
  //! \xmlrestart Individuals are save in \<Individual\> tags in two possible format:
  //!             a long explicit format 
  //! \code
  //    <Individual>
  //      <Structure>
  //        ...
  //      </Structure>
  //    </Individual>
  //! \endcode 
  //!             where the \<Structure\> tags contain the complete information
  //!             to the structure. a compact format:
  //! \code
  //    <Individual string="010011100" />
  //! \endcode 
  //!             where the bitstring objects are encoded as such.
  template< class T_INDIVIDUAL >
  class Evaluator : public GA::Evaluator< T_INDIVIDUAL >
  {
    //! \cond
    public:
      typedef T_INDIVIDUAL t_Individual;

    protected:
      typedef typename t_Individual ::t_IndivTraits t_IndivTraits;
      typedef typename t_IndivTraits::t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_FourierRtoK t_FourierRtoK;
      typedef typename t_IndivTraits::t_FourierKtoR t_FourierKtoR;
      typedef GA::Evaluator<t_Individual> t_Base;
      typedef Evaluator<t_Individual>     t_This;

    protected:
      using t_Base :: current_object;

    public:
      using t_Base::Load;
    //! \endcond

    protected:
      //! The lattice for which decoration search is done
      Ising_CE :: Lattice lattice;
      //! The structure (cell-shape) for which decoration search is done
      Ising_CE :: Structure structure;
      //! The epitaxial growth direction
      atat::rVector3d direction;
      //! The number of independent layers
      types::t_unsigned multiplicity;
      //! \brief A concentration instance as define by
      //!        T_INDIVIDUAL::t_IndivTraits::t_Concentration
      t_Concentration concentration;

    public:
      //! Constructor
      Evaluator() : t_Base() {}
      //! Copy Constructor
      Evaluator   ( const Evaluator &_c )
                : t_Base(), lattice(_c.lattice), structure( _c.structure),
                  direction(_c.direction), multiplicity(_c.multiplicity),
                  concentration(_c.concentration) {}
      //! Destructor
      ~Evaluator() {}

      //! \brief Loads the lattice, the epitaxial parameters from XML, and constructs
      //!        the structure.
      bool Load( const TiXmlElement &_node );
      //! Loads an individual from XML.
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      //! Saves an individual to XML.
      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;

      //! \brief Allows %GA to use physical operators, as define in gaoperator.h
      //! \return a pointer to a %GA operator. This pointer is owned by the callee.
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el )
       { return GA::LoadGaOp<t_Individual>( _el, structure, concentration ); }
      //! Creates a GA::xTaboo instance if requested.
      //! \return a pointer to a Taboo_Base functor. This pointer is owned by
      //!         the callee.
      GA::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );
      //! Creates random individuals using GA::Random.
      bool initialize( t_Individual &_indiv );
      //! Checks from \<GA\> attributes wether to fix the concentration.
      void LoadAttribute ( const TiXmlAttribute &_att )
        { concentration.LoadAttribute( _att ); };

      //! Sets Layered::Evaluator::Structure to \a _indiv
      void init( t_Individual &_indiv );

      //! Prints conceration attributes and structure
      std::string print() const;

    protected:
      //! Loads epitaxial growth parameters and constructs the structure
      bool Load_Structure( const TiXmlElement &_node );
  };

  //! \brief Strict Weak Ordering functor according to depth along eptiaxial
  //!        direction
  //! \details Two vectors are compared by using the value of their scalar
  //           product with Depth::a0. If these scalar product are equal (as
  //           defined by Fuzzy::eq()), then their
  //           scalar product with Depth::a1 are compared. If again these are
  //           equal, then the scalar porducts with Depth::a2 are compared and
  //           the result return. Depth::a0 is the first column of the matrix
  //           given in the constructor argument. Depth::a2 is the second
  //           column, and Depth::a3 the third.
  class Depth
  {
    protected:
      atat::rVector3d a0; //!< First ordering direction
      atat::rVector3d a1; //!< Second ordering direction
      atat::rVector3d a2; //!< Third ordering direction

    public:
      //! \brief Constructor and Initializer
      //! \param _mat Depth::a0 is set to the first column of this matrix,
      //!             Depth::a1 to the second, and Depth::a2 to the third.
      Depth   ( const atat::rMatrix3d &_mat )
            : a0(_mat.get_column(0)), a1(_mat.get_column(1)),
              a2(_mat.get_column(2)) {}
      //! Copy Constructor.
      Depth( const Depth &_c) : a0(_c.a1), a1(_c.a1), a2(_c.a2) {}

      //! Strict weak ordering operator.
      bool operator()(const atat::rVector3d& _1, const atat::rVector3d& _2 );
  };


  //! \brief Taboo functor which forbids large layers. 
  //! \details This is meant for epitaxial growth. Large layers will generally
  //!          accumulate strain unti at some point the surface explodes. The
  //!          point of this taboo is to fordid the appearance of layered
  //!          structure with large homogeneous layers.
  template< class T_INDIVIDUAL >
  class Taboo : public GA::Taboo_Base< T_INDIVIDUAL >
  {
    public:
      //! All relevant %GA %types.
      typedef T_INDIVIDUAL t_Individual;

    protected:
      //! Type of the base class
      typedef GA::Taboo_Base< t_Individual > t_Base;

    protected:
      //! \brief Each component correspond to an atom and is true if the atom
      //!        type is not frozen.
      std::vector<types::t_int> sites;
      //! The maximum number of layers for type "-1" atoms.
      types::t_unsigned d0;
      //! The maximum number of layers for type "1" atoms.
      types::t_unsigned d1;
      //! \brief The largest layer that the current cell-shape can accomodate.
      //! \details This is a geometrical constraint, not a physical constraint.
      types::t_unsigned Nmax;

    public: 
      //! Constructor and Initializer
      Taboo ( const Ising_CE::Structure &_str );

      //! Copy Constructor
      Taboo   ( const Taboo &_c )
            : t_Base( _c ), sites( _c.sites ), d0( _c.d0 ), d1 ( _c.d1 ) {}
      //! Destructor
      ~Taboo() {}

      //! Loads the maximum number of layers from XML.
      bool Load( const TiXmlElement &_node );

      //! Returns true of the individual \a _indiv is taboo. Eg finds a heavy layer.
      bool operator()( const t_Individual &_indiv ) const;
  };


} // namespace Layered
/** @} */

#include "layered.impl.h"

#endif // _LAYERED_H_
