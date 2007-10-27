//
//  Version: $Id$
//
#ifndef _MOLECULARITY_H_
#define _MOLECULARITY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <ostream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <vff/layered.h>
#include <pescan_interface/interface.h>
#include <lamarck/structure.h>
#include <opt/opt_function_base.h>
#include <opt/opt_minimize_gsl.h>
#include <opt/types.h>

#include "layered.h"
#include "pescan.h"
#include "vff.h"
#include "evaluator.h"
#include "individual.h"
#include "gaoperators.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Molecularity
{
  types::t_real inplane_stress( const atat::rMatrix3d &_stress, const atat::rVector3d &_dir );

  class Object : public Layered::Object<>, public Vff::Keeper, public Pescan::Keeper
  {
    protected:
      typedef Layered :: Object<> t_LayeredBase;
      typedef Vff::Keeper         t_VffBase;
      typedef Pescan::Keeper      t_PescanBase;

    public:
      typedef t_LayeredBase :: t_Type t_Type;
      typedef t_LayeredBase :: t_Container t_Container;

    public:
      Object() {}
      Object(const Object &_c) : t_LayeredBase(_c), t_VffBase(_c), t_PescanBase(_c) {};
      ~Object() {};
      
      bool Load( const TiXmlElement &_node )
        { return t_VffBase::Load(_node) and t_PescanBase::Load(_node); }
      bool Save( TiXmlElement &_node ) const
        { return t_VffBase::Save(_node) and t_PescanBase::Save(_node); }
  };

  //! \brief Explicitely defines stream dumping of Molecularity::Object 
  //! \details Modulates the print out for all formats but XML. 
  //! \warning don't touch the dumping to a string, unless you want to change
  //!          read/write to XML. It's got nothing to do with anything here...
  //!          Just repeating myself.
  std::ostream& operator<<(std::ostream &_stream, const Object &_o);
  

  typedef Individual::Types< Object, 
                             Layered::Concentration<2>, 
                             Layered::Fourier<2>    > :: Vector t_Individual;

  class Evaluator : public Layered::Evaluator< t_Individual >
  {
    public:
      typedef Molecularity::t_Individual     t_Individual;
      typedef Traits::GA< Evaluator >        t_GATraits;
    protected:
      typedef Evaluator                      t_This;
      typedef Layered::Evaluator<t_Individual> t_Base;
      typedef Ising_CE::Structure::t_kAtoms  t_kvecs;
      typedef Ising_CE::Structure::t_Atoms   t_rvecs;

    public:
      using t_Base :: Load; 

    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      Pescan::Darwin pescan;
      Vff::Darwin<Vff::Layered> vff;

    public:
      Evaluator() : t_Base(), pescan(structure), vff(structure) {}
      Evaluator   ( const Evaluator &_c )
                : t_Base(_c), pescan(_c.pescan), vff(_c.vff) {}
      ~Evaluator() {}

      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
       { return _indiv.Object().Save(_node) and t_Base::Save( _indiv, _node, _type ); }
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );

      void evaluate();
      eoF<bool>* LoadContinue(const TiXmlElement &_el );

    protected:
      void object_to_quantities( t_Individual & _indiv );
  };

}

#include "molecularity.impl.h"

#endif // _MOLECULARITY_H_
