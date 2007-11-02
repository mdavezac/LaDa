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
#include <opt/opt_function_base.h>
#include <opt/opt_minimize_gsl.h>
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

  //! Object
  struct Object : public TwoSites::Object, public Pescan::Keeper
  {
    friend std::ostream& operator<<(std::ostream &_stream, const Object &_o);
    typedef TwoSites::Object :: t_Container t_Container;
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Object>(BandGap::Object &);
#endif
    types::t_real x, y;


    Object() : TwoSites::Object(), Pescan::Keeper() {}
    Object   (const Object &_c)
           : TwoSites::Object(_c), Pescan::Keeper(_c),
             x(_c.x), y(_c.y) {};
    ~Object() {};
  };

  inline std::ostream& operator<<(std::ostream &_stream, const Object &_o)
  { 
    if( _o.Container().size() <= 30 )
      _stream << (const SingleSite::Object& ) _o << " ";
    _stream << "x=" << _o.x << " y="  << _o.y 
            << (const Pescan::Keeper&) _o;
    return _stream; 
  } 

  typedef Individual::Types< BandGap::Object, 
                             TwoSites::Concentration, 
                             TwoSites::Fourier        > :: Scalar t_Individual;

  class Evaluator : public TwoSites::Evaluator< BandGap::t_Individual >
  {
    public:
      typedef BandGap::t_Individual t_Individual;
      typedef Traits::GA< Evaluator > t_GATraits;
    protected:
      typedef Evaluator t_This;
      typedef TwoSites::Evaluator< t_Individual > t_Base;
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;

    public:
      using t_Base :: Load;
      using t_Base :: Save;
    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      Pescan::Darwin pescan;
      Vff::Darwin<Vff::Functional> vff;

    public:
      Evaluator() : t_Base(), pescan(structure), vff(structure) {}
      Evaluator   ( const Evaluator &_c )
                : t_Base(_c), pescan(_c.pescan), vff(_c.vff) {}
      ~Evaluator() {};

      bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      void LoadAttribute ( const TiXmlAttribute &_att ) {};
      eoF<bool>* LoadContinue(const TiXmlElement &_el )
        { return new GA::mem_zerop_t<Pescan::Darwin>( pescan, &Pescan::Darwin::Continue,
                                                      "Pescan::Continue" );     }

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
