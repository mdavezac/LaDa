//
//  Version: $Id$
//
#ifndef _TWOSITES_H_
#define _TWOSITES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <algorithm>
#include <functional>
#include <string>
#include <sstream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>

#include "evaluator.h"
#include "concentration.h"
#include "functors.h"
#include "taboos.h"
#include "gatraits.h"
#include "single_site.h"
#include "gaoperators.h"

/** \ingroup Genetic
 * @{ */
//! \brief Defines base classes for ternary and quaternary semi-conductors
//! \details It is expected the lattice will contain two different sites, one
//!          for cations, another for anions.  The atoms in the structure are
//!          rearranged (see TwoSites::rearrange_structure) such that first
//!          lattice site of unit cell is immediatly followed by the second
//!          lattice site. The concentration of the lattice sites can be
//!          optimized independantly, or linked through some load-minimizing
//!          function (see X_vs_Y ). The fourier transform assign +/-1 to the
//!          first lattice site and +/- i to the second lattice site.
namespace TwoSites
{
  //! \brief rearranges the atoms of Crystal::Structure such that the fist
  //!        lattice-site is always followed by the second lattice-site.
  void rearrange_structure(Crystal::Structure &);

  //! The type of the object
  typedef SingleSite :: Object Object;

  //! \brief Defines fourier transforms for two-site lattices
  //! \details The fourier transform assign +/-1 to the
  //!          first lattice site and +/- i to the second lattice site.
  struct Fourier
  {
    //! \brief From real to k space
    //! \param[in, out] _rfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type)
    //! \param[in, out] _rend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type< std::complex >)
    template<class T_R_IT, class T_K_IT>
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend );
    //! \brief From k space to real space. 
    //! \details The first range is most likely some instance of
    //!          Crystal::Structure::t_Atoms. The second range is similarly an
    //!          instance of Crystal::Structure::t_kAtoms. The third range
    //!          should be iterators to std::complex.
    //! \pre The range [ \a _rout, \a _rout += \a _rfirst - \a _rend  ) should
    //!      be valid.
    //! \param[in] _rfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type)
    //! \param[in] _rend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type)
    //! \param[in] _kfirst iterator to the first real space atom (of a type
    //!                similar to Crystal::Atom_Type< std::complex >)
    //! \param[in] _kend iterator to the last real space atom (of a type
    //!              similar to Crystal::Atom_Type< std::complex >)
    //! \param[out] _rout iterator to the first complex real-space
    //!              occupation ( of std::complex type )
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend,
             T_O_IT _rout ); // sets rvector values from kspace values
  };

  //! \brief All concentration related behaviors
  //! \details Works for a single, set cell-shape. The concentration of the two
  //!          sites can be fixed. They can also be linked through a
  //!          load-minimizing function (see X_vs_y).  Finally, they can be
  //!          independant. In cases where the concentration is fully or
  //!          partially constrained, this class is able to set the right
  //!          concentration of an Crystal::Structure instance, and of a
  //!          TwoSites::Object instance.  It works correctly with structure
  //!          for which the occupation of some sites are frozen.
  //! \xmlinput see X_vs_y.
  class Concentration : public X_vs_y
  {
    public:
      //! Concentration of the first site. Output variable for Concentration::get()
      types::t_real x;
      //! Concentration of the second site. Output variable for Concentration::get()
      types::t_real y;
      //! The number of sites in this cell-shape.
      types::t_unsigned N;
      //! The number of frozen first lattice sites
      types::t_int Nfreeze_x;
      //! The number of frozen second lattice sites
      types::t_int Nfreeze_y;
      //! \brief Record wether the components of a Object::container is a first or
      //!        second lattice site.
      std::vector<bool> sites;

    public:
      //! Constructor
      Concentration  () 
                    : X_vs_y(), x(0), y(0), N(0),
                      Nfreeze_x(0), Nfreeze_y(0) {}
      //! Copy Constructor
      Concentration   ( const Concentration &_conc)
                    : X_vs_y(_conc), x(_conc.x), y(_conc.y),
                      N(_conc.N), Nfreeze_x(_conc.Nfreeze_x), Nfreeze_y(_conc.Nfreeze_y),
                      sites(_conc.sites) {}
      //! Destructor
      ~Concentration() {}

      //! \brief Loads the required behavior from XML (eg constrained,
      //!        partially constrained...)
      bool Load( const TiXmlElement &_node );

      //! \brief Normalizes the site occupations as given by the \b k-vectors. 
      //! \details A "complex" real-space occupation is computed from the
      //!          k-space intensities. For first-lattice sites, normalized
      //!          site-occupations are set to +/- 1 depending on wether the
      //!          complex real-space occupation is on left or right complex
      //!          half-plane. For second-lattice sites, normalized
      //!          site-occupations are set to +/- 1 depending on wether the
      //!          complex real-space occupation is on upper or lower complex
      //!          half-plane. These behaviors are in-line with the encoding
      //!          expected by the fourier transform.  If the concentration is
      //!          fixed, those sites for which the real value of the complex
      //!          occupation  are closest to zero are flipped first.
      //! \see GA::Krossover, GA::KRandom, GA::KMutation.
      void operator()( Crystal::Structure &_str );
      //! \brief Sets the concentration of \a _obj  if the concentration
      //!        is fully or partially constrained
      void operator()( Object &_obj );
      //! \brief Computes the concentration of \a _str and stores the result in
      //!        Concentration::x and Concentration::y
      void get( const Crystal::Structure &_str);
      //! \brief Computes the concentration of \a _obj and stores the result in
      //!        Concentration::x and Concentration::y
      void get( const Object &_obj );
      //! \brief Computes the number of sites in \a _str for which the
      //!        occupation is fixed
      void setfrozen ( const Crystal::Structure &_str );

    protected:
      //! \brief Actually does the job of normalizing the site occupations for
      //!        Concentration::operator()(Crystal::Structure&)
      void normalize( Crystal::Structure &_str, const types::t_int _site, 
                      types::t_real _tochange);

  };

  //! \brief Partially overrides GA::Evaluator default to implement behaviors
  //!        relevant to two-site lattices with a fixed cell-shape.
  //! \details Mostly allows for physical %GA operators (GA::Krossover,
  //!          GA::xTaboo, GA::Crossover, etc...) and links them with the
  //!          bitstring objects defined above, as well as the concentration
  //!          functor TwoSites::Concentration. Saving and Restarting of
  //!          individuals is partially implementated as relates to
  //!          BitString::Object and Crystal::Structure.
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
  template<class T_INDIVIDUAL >
  class Evaluator : public GA::Evaluator< T_INDIVIDUAL >
  {
    //! \cond
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef Traits::GA< Evaluator<t_Individual> > t_GATraits;
    protected:
      typedef typename t_Individual::t_IndivTraits t_IndivTraits;
      typedef typename t_IndivTraits::t_Object t_Object;
      typedef typename t_IndivTraits::t_Concentration t_Concentration;
      typedef GA::Evaluator<t_Individual> t_Base;
      typedef Evaluator<t_Individual> t_This;

    public:
      using t_Base :: Load;

    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      typedef Crystal::Structure::t_kAtoms t_kvecs;
      typedef Crystal::Structure::t_Atoms t_rvecs;
    //! \endcond
      
    protected:
      //! The lattice for which decoration search is done
      Crystal::Lattice lattice;
      //! The structure (cell-shape) for which decoration search is done
      Crystal::Structure structure;
      //! \brief A concentration instance as define by
      //!        T_INDIVIDUAL::t_IndivTraits::t_Concentration
      t_Concentration concentration;

    public:
      //! Constructor
      Evaluator   () {}
      //! Copy Constructor
      Evaluator   ( const t_This &_c )
                : t_Base( _c ),
                  lattice( _c.lattice ),
                  structure( _c.structure ),
                  concentration( _c.concentration ) {}
      //! Destructor
      ~Evaluator() {};


      //! Saves an individual to XML.
      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      //! Loads an individual from XML.
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      //! \brief Loads the lattice and the structure (cell-shape) from XML.
      bool Load( const TiXmlElement &_node );
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
      //! Sets TwoSites::Evaluator::Structure to \a _indiv
      void init( t_Individual &_indiv );

      //! \brief Used to submits individuals to history, etc, prior to starting %GA
      //! \details initializes the endopoints of a convex-hull, for instance.
      //! Presubmitted individuals are not put into the population.
      //! \see GA::Evaluator::presubmit()
      void presubmit( std::list<t_Individual> &_pop );

    protected:
      //! Makes input of structure and lattice are coherent.
      bool consistency_check();
  };

} // namespace TwoSites

namespace GA
{
  //! Policies for keeping track of functional evaluations in objects.
  namespace Keepers
  {
    struct ConcTwo
    {
      //! The concentration of the first site.
      types::t_real x;
      //! The concentration of the second site.
      types::t_real y;
      //! Constructor.
      ConcTwo() : x(-2), y(-2) {}
      //! Copy Constructor.
      ConcTwo( const ConcTwo& _c ) : x(_c.x), y(_c.y) {}
      //! Loads from attributes of \a _node.
      bool Load( const TiXmlElement &_node );
      //! Saves as attributes of \a _node.
      bool Save( TiXmlElement &_node ) const;
      //! Serializes concentration class.
      template<class Archive>
        void serialize(Archive & _ar, const unsigned int _version)
          { _ar & x; _ar & y; }
    };
  }
}

#include "two_sites.impl.h"

/* @} */

#endif // _TWOSITES_OBJECT_H_
