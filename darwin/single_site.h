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

#include "lamarck/structure.h"
#include "opt/types.h"

#include "evaluator.h"
#include "functors.h"
#include "taboos.h"
#include "bitstring.h"
#include "gaoperators.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif


// Defines base classes for two site structures

namespace SingleSite
{
  using Ising_CE::Fourier;

  struct Object : public BitString::Object<> 
  {
    // Typedefs
    protected:
      typedef BitString :: Object<> t_Base;
    public:
      typedef t_Base :: t_Type t_Type;
      typedef t_Base :: t_Container t_Container;


    // Member Functions
    public:
      Object() {}
      Object(const Object &_c) : t_Base(_c) {};
      ~Object() {};
    
  };

  void operator<<(Ising_CE::Structure &_str, const Object &_o);
  void operator<<(Object &_o, const Ising_CE::Structure &_str);

  class Concentration 
  {
    protected:
      types::t_real x0;
    public:
      types::t_real x;
      types::t_unsigned N;
      types::t_int Nfreeze;
      bool single_c;

    public:
      Concentration () : x0(0), x(0), N(0), Nfreeze(0), single_c(false) {}
      Concentration   ( const Concentration &_conc)
                    : x0(_conc.x0), x(_conc.x), N(_conc.N), Nfreeze(_conc.Nfreeze) {}
      ~Concentration() {}

      bool Load( const TiXmlElement &_node );
      void LoadAttribute ( const TiXmlAttribute &_att );

      void operator()( Ising_CE::Structure &_str );
      void operator()( Object &_obj );
      void operator()( const Ising_CE::Structure &_str, Object &_object,
                       types::t_int _concx, types::t_int _concy );
      void set( const Ising_CE::Structure &_str);
      void set( const Object &_obj );
      void setfrozen ( const Ising_CE::Structure &_str );

      std::string print() const;

    protected:
      void normalize( Ising_CE::Structure &_str, 
                      types::t_real _tochange);

  };

  template<class T_INDIVIDUAL>
  class Evaluator : public GA::Evaluator< T_INDIVIDUAL >
  {
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
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;
      
    protected:
      Ising_CE::Lattice lattice;
      Ising_CE::Structure structure;
      t_Concentration concentration;

    public:
      Evaluator   ()
                : lattice(), structure(), concentration() {}
      Evaluator   ( const t_This &_c ) 
                : lattice( _c.lattice), structure(_c.structure),
                  concentration(_c.concentration) {}
      ~Evaluator() {};


      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el )
       { return GA::LoadGaOp<t_Individual>( _el, structure, concentration ); }
      GA::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el )
      {
        if ( concentration.single_c ) return NULL;
        GA::xTaboo<t_Individual> *xtaboo = new GA::xTaboo< t_Individual >( concentration );
        if ( xtaboo and xtaboo->Load( _el ) )  return xtaboo;
        if ( xtaboo ) delete xtaboo;
        return NULL;
      }

      bool initialize( t_Individual &_indiv )
      {
        GA::Random< t_Individual > random( concentration, structure, _indiv );
        _indiv.invalidate(); return true;
      }
      void LoadAttribute ( const TiXmlAttribute &_att )
        { concentration.LoadAttribute( _att ); };

      //! \brief Used to submits individuals to history, etc, prior to starting %GA
      //! \details initializes the endopoints of a convex-hull, for instance.
      //! Presubmitted individuals are not put into the population.
      //! \see GA::Evaluator::presubmit()
      void presubmit( std::list<t_Individual> &_pop );


    protected:
      bool consistency_check();
  };

} // namespace TwoSites

#include "single_site.impl.h"

#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize< SingleSite::Object >
                                       ( SingleSite::Object & _object )
  {
    return serialize( _object.bitstring );
  }
}
#endif

#endif // _TWOSITES_OBJECT_H_
