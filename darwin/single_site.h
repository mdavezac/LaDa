//
//  Version: $Id$
//
#ifndef _SINGLE_SITE_H_
#define _SINGLE_SITE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>

#include "evaluator.h"
#include "functors.h"
#include "taboos.h"
#include "bitstring.h"
#include "gaoperators.h"
#include "scaling.h"


/** \ingroup Genetic
 * @{ */
//! \brief Defines \e physical base classes for single-site, single-cell-shape, decoration search.
//! \details By single-site, we mean one site in lattice unit-cell. By
//!          single-cell-shape, we mean that the whole process takes place for
//!          one type of cell-shape (ICS, see <A
//!          HREF="http://dx.doi.org/10.1103/PhysRevB.74.014204" > Trimarchi \e
//!          et \e al., PRB \b 74, p14204, (2006) </A>). By decoration, we mean
//!          that only the lattice occupations are searched (and not lattice
//!          type or others). This is the type of search performed <A
//!          HREF="http://dx.doi.org/10.1088/0953-8984/19/40/402201"> here </A>
//!          with the objective being a simple minimization.
//!          More specifically, the fourier transform from Crystal is used. A
//!          Bitstring object is defined which can load and be loaded from an
//!          Crystal::Structure instance. A concentration class is created
//!          which can ecaluate the contration of SingleSite::Object and of an
//!          Crystal::Structure, whatever the number and type of frozen atoms.
//!          It can also set the concentration of SingleSite::Object and 
//!          Crystal::Structure instances. Finally, and evaluator class is
//!          derived from GA::Evaluator to use physical %GA operators.
//! \xmlinput If the \<GA\> tag contains an x0 attribute, then the
//!           concentration is fixed. The structure which the individuals
//!           characterize is expected to be found directly under the overall
//!           \<Job\> tag.  Its exact format is described in \ref TagStructure.
namespace SingleSite
{
  //! \brief Type of the Fourier functor.
  //! \details Once doxygen is smart enough, should become
  //! \code
  //  using Crystal::Fourier;
  //! \endcode
  //!      
  typedef Crystal::Fourier Fourier;

  //! \brief Object describing the decoration of a structure.
  //! \details Mostly, this is a redefinition of a BitString::Object<>. The
  //!          point is to be able to define 
  //!          SingleSite::operator<<( Crystal::Structure&, const Object&) and 
  //!          SingleSite::operator<<( Object&, const Crystal::Structure&) without
  //!          interfering with possible previous
  //!          declarations.
  //! \see SingleSite::operator<<( Crystal::Structure&, const Object& ), 
  //!      SingleSite::operator<<( Object&, const Crystal::Structure& ). 
  struct Object : public BitString::Object<> 
  {
    protected:
      //! \cond
      typedef BitString :: Object<> t_Base;
      //! \endcond
    public:
      //! \brief See function::Base::t_Type
      typedef t_Base :: t_Type t_Type;
      //! See function::Base::t_Container
      typedef t_Base :: t_Container t_Container;


    public:
      //! Constructor
      Object() {}
      //! Copy Constructor
      Object(const Object &_c) : t_Base(_c) {};
      //! Destructor
      ~Object() {};
    
      //! Serializes a single-site object.
      template<class Archive> void serialize(Archive & _ar, const unsigned int _version)
        { _ar & boost::serialization::base_object< BitString::Object<> >(*this); }
  };

  //! Dumps the decoration of \a _str into the object \a _o
  void operator<<(Crystal::Structure &_str, const Object &_o);
  //! Dumps the decoration \a _o into the structure \a _str
  void operator<<(Object &_o, const Crystal::Structure &_str);

  //! \brief Concentration related behaviors for the decorations of single-site
  //!        structures.
  //! \details It can compute and set the concentrations of Object and
  //!          Crystal::Structure instances for single-site lattices. It works
  //!          correctly with structure for which the occupation of some sites
  //!          are frozen.
  //! \xmlinput Specifying an x0 attribute in the \<%GA\> tag implies that the
  //!           concetration will be set to that value throughout.
  class Concentration 
  {
    protected:
      //! \brief Value of the fixed concentration.
      //! \details Not used if the concentration is not fized
      types::t_real x0;
    public:
      //! \brief Stores the concentation of the last computed Object of
      //!        Crystal::Structure instance..
      types::t_real x;
      //! The number of sites in the current cell-shape
      types::t_unsigned N;
      //! \brief The number of sites \e with \e fixed \e occupations in the
      //!        current cell-shape
      types::t_int Nfreeze;
      //! True if the concentration is fixed
      bool single_c;

    public:
      //! Constructor
      Concentration () : x0(0), x(0), N(0), Nfreeze(0), single_c(false) {}
      //! Copy Constructor
      Concentration   ( const Concentration &_conc)
                    : x0(_conc.x0), x(_conc.x), N(_conc.N), Nfreeze(_conc.Nfreeze) {}
      //! Destructor
      ~Concentration() {}

      //! \brief Loads the \e fixed concentration from a \<Concentration\> tag. 
      //! \details If it is not present, then the concentration is not fixed.
      bool Load( const TiXmlElement &_node );
      //! \brief Chechs whether an attribute is x0. If it is, the concentration
      //!        is fixed.
      //! \details This is used mainly to load attributes from the \<GA\> Tag
      //!          by GA::Darwin.
      void LoadAttribute ( const TiXmlAttribute &_att );

      //! \brief Normalizes the site occupations as given by the \b k-vectors. 
      //! \details A "complex" real-space occupation is computed from the
      //!          k-space intensities. Normalized site-occupations are set to
      //!          +/- 1 depending on which half=plane the complex value are.
      //!          If the concentration is fixed, those sites for which the
      //!          real value of the complex occupation  are closest to zero
      //!          are flipped first.
      //! \see GA::Krossover, GA::KRandom, GA::KMutation.
      void operator()( Crystal::Structure &_str );
      //! \brief Sets \a _obj to concentration Concentration::x0 if the concentration
      //!        is fixed.
      void operator()( Object &_obj );
      //! \brief Computes the concentration of \a _str and stores the result in
      //!        Concentration::x
      void get( const Crystal::Structure &_str);
      //! \brief Computes the concentration of \a _obj and stores the result in
      //!        Concentration::x
      void get( const Object &_obj );
      //! \brief Computes the number of sites in \a _str for which the
      //!        occupation is fixed
      void setfrozen ( const Crystal::Structure &_str );

      //! \brief Returns a string describing the statuc of this instance
      //! \returns "Single Concentration, x0 =?" with ?=Concentration::x0 or 
      //!          "Concentration Range"
      std::string print() const;

    protected:
      //! \brief Actually does the job of normalizing the site occupations for
      //!        Concentration::operator()(Crystal::Structure&)
      void normalize( Crystal::Structure &_str, 
                      types::t_real _tochange);

  };

  //! \brief Partially overrides GA::Evaluator default to implement behaviors
  //!        relevant to single-lattice-site structure decoration search. 
  //! \details Mostly allows for physical %GA operators (GA::Krossover,
  //!          GA::xTaboo, GA::Crossover, etc...) and links them with the
  //!          bitstring objects defined above, as well as the concetration
  //!          functor SingleSite::Concentration. Saving and Restarting of
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
  template<class T_INDIVIDUAL>
  class Evaluator : public GA::Evaluator< T_INDIVIDUAL >
  {
    //! \cond
    public:
      typedef T_INDIVIDUAL t_Individual;
    protected:
      typedef typename t_Individual::t_IndivTraits t_IndivTraits;
      typedef typename t_IndivTraits::t_Object t_Object;
      typedef typename t_IndivTraits::t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_FourierRtoK t_FourierRtoK;
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
      Evaluator   ()
                : lattice(), structure(), concentration() {}
      //! Copy Constructor
      Evaluator   ( const t_This &_c ) 
                : lattice( _c.lattice), structure(_c.structure),
                  concentration(_c.concentration) {}
      //! Destructor
      ~Evaluator() {};


      //! Saves a T_INDIVIDUAL instance to XML
      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      //! Load a T_INDIVIDUAL instance from XML
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      //! \brief Allows %GA to use physical operators, as define in gaoperator.h
      //! \return a pointer to a %GA operator. This pointer is owned by the callee.
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el )
       { return GA::LoadGaOp<t_Individual>( _el, structure, concentration ); }
      //! Creates a GA::xTaboo instance if requested.
      //! \return a pointer to a Taboo_Base functor. This pointer is owned by the callee.
      GA::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el )
      {
        if ( concentration.single_c ) return NULL;
        GA::xTaboo<t_Individual> *xtaboo = new GA::xTaboo< t_Individual >( concentration );
        if ( xtaboo and xtaboo->Load( _el ) )  return xtaboo;
        if ( xtaboo ) delete xtaboo;
        return NULL;
      }

      //! Creates random individuals using GA::Random.
      bool initialize( t_Individual &_indiv )
      {
        GA::Random< t_Individual > random( concentration, structure, _indiv );
        _indiv.invalidate(); return true;
      }
      //! Checks from \<GA\> attributes wether to fix the concentration.
      void LoadAttribute ( const TiXmlAttribute &_att )
        { concentration.LoadAttribute( _att ); };

      //! \brief Used to submits individuals to history, etc, prior to starting %GA
      //! \details initializes the endopoints of a convex-hull, for instance.
      //! Presubmitted individuals are not put into the population.
      //! \see GA::Evaluator::presubmit()
      void presubmit( std::list<t_Individual> &_pop );

      //! \brief Returns a functor for printing xyz format
      //! \see Evaluator::LoadPrintBest()
      eoMonOp<const t_Individual>* LoadPrintBest( const TiXmlElement &_node );

      //! Prints conceration attributes and structure
      std::string print() const;

    protected:
      //! \brief Does simple consistency checks between the structure as loaded from
      //!        \<Structure\> and the lattice, as loaded from \<Lattice\>
      //! \details If the structure is loaded after the lattice, and if
      //!          Crystal::Structure::Lattice is set corretly, this should
      //!          probably not be necessary.
      bool consistency_check();
  };

  /** \brief Defines a phenotypic distance which takes into account the quantities
   *   and the concentration of an individuals
       \f$ \alpha_x | x_{\sigma_i} - x_{\sigma_j} | 
           + \sum_{e=0}^{D}\alpha_e |q_{\sigma_i} - q_{\sigma_j} |  \f$ */
  template<class T_GATRAITS, types::t_unsigned _D = 1>
  class Distance
  {
    public:
      typedef T_GATRAITS t_GATraits; //!< All %GA %types
      //! Stores how many quantities to look at
      const static types::t_int d = _D; 

    protected:
      //! Type of the individuals
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! Type of the object characterising the individuals
      typedef typename t_GATraits :: t_Object t_Object;
      //! Type of the \b scalar fitness
      typedef typename t_GATraits :: t_ScalarFitness t_ScalarFitness;
      //! Type of the  quantity of the \b scalar fitness
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;
      //! All types of the \e physical individual
      typedef typename t_GATraits :: t_IndivTraits  t_IndivTraits;
      //! Type of the concentration functor
      typedef typename t_IndivTraits :: t_Concentration t_Concentration;

    protected:
      //! \brief Reference to the concentration for this individual
      //! \details This is the functor which computes the concentration of an individual
      t_Concentration &concentration;
      //! \brief \f$\alpha_x\f$, factors the concentration
      t_ScalarFitnessQuantity xcoef;
      //! \brief \f$\alpha_e\f$. The _D factors of the first \a _D quantities.
      t_ScalarFitnessQuantity qcoefs[_D];
  
    public:
      //! Constructor
      Distance  ( t_Concentration &_conce ) 
               : concentration( _conce ), xcoef(0)
        { for( unsigned u=0; u < _D; ++u ) qcoefs[u] = 0; }
      //! Copy Constructor
      Distance  ( const Distance &_c ) 
               : concentration( _c.concentration ), xcoef(_c.xcoef)
        { for( unsigned u=0; u < _D; ++u ) qcoefs[u] = _c.qcoefs[u]; }
      /** Returns \f$ \alpha_x | x_{\sigma_i} - x_{\sigma_j} | 
                      + \sum_{e=0}^{D}\alpha_e |q_{\sigma_i} - q_{\sigma_j} |  \f$ */
      t_ScalarFitnessQuantity operator()( const t_Individual &_i1,
                                          const t_Individual &_i2) const;
      //! Returns "SingleSite::Distance"
      std::string what_is() const;

      //! Loads the phenotypic distance
      bool Load( const TiXmlElement &_node );
  }; 

  //! Loads a Niche< Triangular::Sharing< Distance > > from XML
  template<class T_GATRAITS, types::t_int _D> Scaling::Base<T_GATRAITS>*
    new_Niche_from_xml( const TiXmlElement &_node,
                        typename T_GATRAITS :: t_IndivTraits :: t_Concentration &_conce );

  template<class T_INDIVIDUAL> 
    inline std::string Evaluator<T_INDIVIDUAL> :: print() const
    {
      std::ostringstream sstr;
      sstr << concentration.print() << "\n" << structure;
      return sstr.str();
    }
} // namespace TwoSites
/* @} */

#include "single_site.impl.h"

#endif // _SINGLESITES_H_
