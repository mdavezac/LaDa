//
//  Version: $Id$
//
#ifndef _BANDGAP_H_
#define _BANDGAP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <vff/functional.h>
#include <pescan_interface/interface.h>
#include <lamarck/structure.h>
#include <opt/function_base.h>
#include <opt/gsl_minimizers.h>
#include <opt/types.h>

#include "two_sites.h"
#include "evaluator.h"
#include "concentration.h"
#include "functors.h"
#include "individual.h"
#include "vff.h"
#include "pescan.h"


#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

/** \ingroup Genetic
 * @{ */
//! \brief Implements single-cell-shape decoration search for band gaps.
//! \details The strain is relaxed using the Valence Force Field developped in
//!          namespace Vff. The band-gap is obtained from pescan \e via the
//!          Pescan::Interface class. 
//! 
//!          The band gap can be obtained from a complete all-electron
//!          diagonalization or, more efficiently, through a partial
//!          diagonalization using the folded spectrum method (see <A
//!          HREF="http://dx.doi.or/10.1103/PhysRevB.51.17398"> L-W Wang and A.
//!          Zunger PRB \b 51, 17398 (1995) </A>). The latter needs two
//!          reference energies to do find the band-gap. If these are \b not
//!          given on input, Pescan::Darwin will first perform a full
//!          diagonalization and peg the reference energies to the HOMO and the
//!          LUMO. From there, only Folded-Spectrum methods are performed
//!          unless one of two situations arise: 
//!            - A metallic "band gap" is found 
//!            - The user has requested all-electron calculations every \a N
//!              generations
//!            .
//! \see For more details on the functionals, see Pescan::Interface and
//!      Vff::Functional
namespace BandGap
{

  //! \brief BitString Object with Pescan capacity
  //! \details Other than the bitstring itself and the Pescan::Keeper
  //           variables, this object also stores the x and y concentrations of
  //           a quaternary. It overloads dumping an object to a stream.
  //! \see BandGap::operator<<( std::ostream &, const Object& )
  struct Object : public TwoSites::Object, public Pescan::Keeper
  {
    //! The type of the BitString container
    typedef TwoSites::Object :: t_Container t_Container;
    //! The concentration of the first site.
    types::t_real x;
    //! The concentration of the second site.
    types::t_real y;

    //! Constructor
    Object() : TwoSites::Object(), Pescan::Keeper() {}
    //! Copy Constructor
    Object   (const Object &_c)
           : TwoSites::Object(_c), Pescan::Keeper(_c),
             x(_c.x), y(_c.y) {};
    //! Destructor
    ~Object() {};
  };

  //! \brief Dumps a BandGap::Object to a stream.
  //! \details  This routine is used for printing results only, and never to
  //!           serialize in XML format.
  inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  { 
    if( _o.Container().size() <= 30 )
      _stream << (const SingleSite::Object& ) _o << " ";
    _stream << "x=" << _o.x << " y="  << _o.y 
            << (const Pescan::Keeper&) _o;
    return _stream; 
  } 

  //! \brief Type of the \e physical BandGap individual
  //! \details In addition to BandGap object, the individual uses
  //!          TwoSites::Concentration and TwoSites Fourier. By default, it
  //!          declares a \e scalar fitness.
  typedef Individual::Types< BandGap::Object, 
                             TwoSites::Concentration, 
                             TwoSites::Fourier > :: Scalar t_Individual;

  //! \brief %Evaluator class for band-gap decoration search
  //! \details A Pescan::Darwin and a Vff::Darwin<Vff::Functional> objects are
  //!          declared which allow for strain and bandgap computation.
  //!          The \e physical ga operators of gaoperators.h are used to mate
  //!          individuals. \e A \e priori, I can't see what we would need to
  //!          compute a convex-hull for, so GA::Evaluator::presubmit() is
  //!          overriden to simply clear the population in its argument.
  class Evaluator : public TwoSites::Evaluator< BandGap::t_Individual >
  {
    public:
      //! \brief The type of the \e physical individual used in this decoration
      //! search.
      typedef BandGap::t_Individual t_Individual;
      //! All pertinent %GA traits
      typedef Traits::GA< Evaluator > t_GATraits;
    protected:
      //! \cond
      typedef Evaluator t_This;
      typedef TwoSites::Evaluator< t_Individual > t_Base;
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;
      //! \endcond

    public:
      using t_Base :: Load;
      using t_Base :: Save;
    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      //! interface to... uh ... the pescan interface.
      Pescan::Darwin pescan; 
      //! Interface to Vff::Functional.
      Vff::Darwin<Vff::Functional> vff; 

    public:
      //! Constructor
      Evaluator() : t_Base(), pescan(structure), vff(structure) {}
      //! Copy Constructor
      Evaluator   ( const Evaluator &_c )
                : t_Base(_c), pescan(_c.pescan), vff(_c.vff) {}
      //! Destructor
      ~Evaluator() {};

      //! Saves an individual to XML
      bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      //! Loads an individual from XML
      bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      //! Loads structure, lattice, pescan, vff from XML
      bool Load( const TiXmlElement &_node );
      //! Allows pescan all-electron recomputation.
      eoF<bool>* LoadContinue(const TiXmlElement &_el )
        { return new GA::mem_zerop_t<Pescan::Darwin>( pescan, &Pescan::Darwin::Continue,
                                                      "Pescan::Continue" );     }

      //! Evaluates the band gap after strain minimization
      void evaluate();

      //! \brief Pescan is costly and I'm not sure why we would want to do a
      //!        convex-hull, so presubmit does nothing in this implementation.
      //! Presubmitted individuals are not put into the population.
      //! \see GA::Evaluator::presubmit(), TwoSites::Evaluator::presubmit()
      void presubmit( std::list<t_Individual> &_pop ) { _pop.clear(); }
  };


} // namespace BandGap

/** @} */

#ifdef _MPI
namespace mpi
{
  /** \ingroup Genetic
   *  \brief Serializes a BandGap::Object.
   *  \details Includes the serialization of the BitString::container, the
   *  concentrations x and y, and the Pescan::Keeper member variables. */
  template<>
  inline bool mpi::BroadCast::serialize<BandGap::Object>( BandGap::Object & _object )
  {
    return     serialize<Pescan::Keeper>( _object )
           and serialize( _object.x )
           and serialize( _object.y )
           and _object.broadcast( *this );
  }
}
#endif

#endif // _BANDGAP_H_
