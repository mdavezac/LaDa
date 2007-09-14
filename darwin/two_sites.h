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

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include "lamarck/structure.h"
#include "opt/types.h"

#include "evaluator.h"
#include "concentration.h"
#include "functors.h"
#include "taboos.h"
#include "gatraits.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "single_site.h"


// Defines base classes for two site structures

namespace TwoSites
{
  void rearrange_structure(Ising_CE::Structure &);

  typedef SingleSite :: Object Object;

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
  class Evaluator : public darwin::Evaluator< T_INDIVIDUAL, T_INDIV_TRAITS >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIV_TRAITS t_IndivTraits;
    protected:
      typedef typename t_IndivTraits::t_Object t_Object;
      typedef darwin::Evaluator<t_Individual, t_IndivTraits> t_Base;
      typedef Evaluator<t_Individual, t_IndivTraits> t_This;

    public:
      using t_Base :: Load;
      using t_Base :: Save;
    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;
      
    protected:
      Ising_CE::Lattice lattice;
      Ising_CE::Structure structure;
      types::t_real crossover_probability;
      types::t_real x, y;
      types::t_real lessthan, morethan;
      X_vs_y x_vs_y;

    public:
      Evaluator   ()
                : crossover_probability(0.5), 
                  x(0), y(0), lessthan(1.0), morethan(-1.0) {}
      Evaluator   ( const t_Base &_c )
                : crossover_probability( &_c.crossover_probability ), 
                  x(_c.x), y(_c.y), lessthan(_c.lessthan), morethan(_c.moerethan) {}
      ~Evaluator() {};


      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el );
      darwin::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );

      bool Krossover( t_Individual  &_offspring, const t_Individual &_parent,
                      bool _range = false );
      bool Crossover( t_Individual &_indiv1, const t_Individual &_indiv2 );

      bool initialize( t_Individual &_indiv );

    protected:
      void get_xy_concentrations( const Ising_CE::Structure &_str );
      bool consistency_check();
      void set_concentration( Ising_CE::Structure &_str );
      void normalize( Ising_CE::Structure &_str, 
                      const types::t_int _site, types::t_real _tochange);
      bool Taboo(const t_Individual &_indiv );
  };

} // namespace TwoSites

#include "two_sites.impl.h"

#endif // _TWOSITES_OBJECT_H_
