//
//  Version: $Id$
//
#ifndef _DARWIN_EMASS_H_
#define _DARWIN_EMASS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <vff/va.h>
#include <lamarck/structure.h>
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

namespace eMassSL
{

  //! \brief BitString Object with effective mass capacity
  //! \details Other than the bitstring itself and the eMassSL::Keeper
  //           variables, this object also stores the x and y concentrations of
  //           a quaternary. It overloads dumping an object to a stream.
  //! \see eMassSL::operator<<( std::ostream &, const Object& )
  struct Object : public Layered::Object<>, 
                  public Vff::Keeper,
                  public eMassSL::Keeper
  {
    protected:
      //! Layered bitstring base
      typedef Layered :: Object<> t_LayeredBase;  
      //! Strain info base 
      typedef Vff::Keeper         t_VffBase;      
      //! Band gap info base
      typedef Pescan::Keeper      t_PescanBase;   

    public:
      //! see function::Base::t_Type
      typedef t_LayeredBase :: t_Type t_Type; 
      //! see function::Base::t_Container
      typedef t_LayeredBase :: t_Container t_Container;

    public:
      //! Constructor
      Object() : t_LayeredBase(), t_VffBase(), t_PescanBase() {}
      //! Copy Constructor
      Object(const Object &_c) : t_LayeredBase(_c), t_VffBase(_c), t_PescanBase(_c) {};
      //! Destructor
      ~Object() {};
      
      //! Loads strain and band-gap info from XML
      bool Load( const TiXmlElement &_node )
        { return t_VffBase::Load(_node) and t_PescanBase::Load(_node); }
      //! Saves strain and band-gap info to XML
      bool Save( TiXmlElement &_node ) const
        { return t_VffBase::Save(_node) and t_PescanBase::Save(_node); }
  };

  //! \brief Explicitely defines stream dumping of eMassSL::Object 
  //! \details Modulates the print out for all formats but XML. 
  //! \warning don't touch the dumping to a string, unless you want to change
  //!          read/write to XML. It's got nothing to do with anything here...
  //!          Just repeating myself.
  std::ostream& operator<<(std::ostream &_stream, const Object &_o);
  
  //! \brief %Individual type for eMassSL.
  //! \details The object type is the one above, eg a BitString::Object adapted
  //!          for Layered structures and containing info for stress and
  //!          band-gap. The concentration functor, as well as the Fourier
  //!          transform functors are also specialized for layered objects.
  typedef Individual::Types< Object, 
                             Layered::Concentration<2>, 
                             Layered::Fourier<2>    > :: Vector t_Individual;


  //! \brief Evaluator class for effective mass search of a layered structure.
  //! \details Mostly, this class defines  a eMassSL::Darwin instance, and a
  //!          Vff::Darwin<Vff::Layered> instance for evaluating (and
  //!          minimizing) in-plane-stress and for evaluating band-gaps.
  class Evaluator : public Layered::Evaluator< t_Individual >
  {
    public:
      //! Type of the individual
      typedef Molecularity::t_Individual     t_Individual;
      //! All %types relevant to %GA
      typedef Traits::GA< Evaluator >        t_GATraits;
    protected:
      //! Type of this class
      typedef Evaluator                      t_This;
      //! Type of the base class
      typedef Layered::Evaluator<t_Individual> t_Base;

    public:
      using t_Base :: Load; 

    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      //! The pescan interface object for obtaining band-gaps
      eMassSL::Darwin pescan;
      //! The vff object for minimizing and computing strain/stress.
      Vff::Darwin<Vff::Layered> vff;

    public:
      //! Constructor
      Evaluator() : t_Base(), pescan(structure), vff(structure) {}
      //! Copy Constructor
      Evaluator   ( const Evaluator &_c )
                : t_Base(_c), pescan(_c.pescan), vff(_c.vff) {}
      //! Destructor
      ~Evaluator() {}

      //! Saves an individual to XML
      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
       { return _indiv.Object().Save(_node) and t_Base::Save( _indiv, _node, _type ); }
      //! Load an individual from XML
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      //! Loads the lattice, layered structure, pescan, and vff from XML.
      bool Load( const TiXmlElement &_node );

      //! Computes the band-gap and in-plane-stress of the current_individual.
      void evaluate();
      //! Allows for periodic all-electron computations
      eoF<bool>* LoadContinue(const TiXmlElement &_el );

    protected:
      //! Transforms stress and band-edges to quantities in \a _indiv.
      void object_to_quantities( t_Individual & _indiv );
  };


} // namespace BandGap

#include "molecularity.impl.h"
/** @} */
#endif // _BANDGAP_H_
