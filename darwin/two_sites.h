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

  class Concentration : public X_vs_y
  {
    public:
      types::t_real x, y;

    public:
      Concentration() : X_vs_y(), x(0), y(0) {}
      Concentration( const Concentration &_conc) : X_vs_y(_conc), x(_conc.x), y(_conc.y) {}
      ~Concentration() {}

      bool Load( const TiXmlElement &_node )
      {
        if( not X_vs_y::Load( _node ) ) return false;
        if( not is_singlec() )  return true;
        x = get_x();  y = get_y();
        return true;
      }

      void operator()( Ising_CE::Structure &_str );
      void operator()( const Ising_CE::Structure &_str, Object &_object,
                       types::t_int _concx, types::t_int _concy );
      void set( const Ising_CE::Structure &_str);

    protected:
      void normalize( Ising_CE::Structure &_str, const types::t_int _site, 
                      types::t_real _tochange);

  };

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
      types::t_real lessthan, morethan;
      Concentration concentration;

    public:
      Evaluator   ()
                : lessthan(1.0), morethan(-1.0) {}
      Evaluator   ( const t_Base &_c )
                : lessthan(_c.lessthan), morethan(_c.moerethan) {}
      ~Evaluator() {};


      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el );
      darwin::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );

      bool initialize( t_Individual &_indiv );

    protected:
      void get_xy_concentrations( const Ising_CE::Structure &_str );
      bool consistency_check();
      void set_concentration( Ising_CE::Structure &_str );
      void normalize( Ising_CE::Structure &_str, 
                      const types::t_int _site, types::t_real _tochange);
      bool Taboo(const t_Individual &_indiv );
  };

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
  class Krossover : public eoGenOp<T_INDIVIDUAL> 
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIV_TRAITS t_IndivTraits;
    protected:
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      Concentration &concentration;
      Ising_CE::Structure &structure;
      types::t_real rate;
      bool do_range;

    public:
      Krossover   ( Concentration &_c, Ising_CE::Structure &_str )
                : concentration(_c), structure(_str), rate(0.5), do_range(false) {}
      Krossover   ( const Krossover &_k )
                : concentration(_k.concentration), structure(_k.structure),
                  rate(_k.rate), do_range(_k.do_range) {}
      ~Krossover() {}

      bool Load( const TiXmlElement &_node );

      virtual std::string className() const { return "TwoSites::Krossover"; }
      unsigned max_production(void) { return 1; } 

      void apply(eoPopulator<t_Individual>& _pop);
      std::string print_out() const
      {
        std::ostringstream sstr;
        sstr << "Krossover rate = " << rate;
        if (do_range) sstr << ", Range = true ";
        return sstr.str();
      }
  };

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
  class KMutation : public eoGenOp<T_INDIVIDUAL> 
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIV_TRAITS t_IndivTraits;
    protected:
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      Concentration &concentration;
      Ising_CE::Structure &structure;
      types::t_real rate;

    public:
      KMutation   ( Concentration &_c, Ising_CE::Structure &_str )
                : concentration(_c), structure(_str), rate(0.5) {}
      KMutation   ( const KMutation &_k )
                : concentration(_k.concentration), structure(_k.structure),
                  rate(_k.rate) {}
      ~KMutation() {}

      bool Load( const TiXmlElement &_node );

      virtual std::string className() const { return "TwoSites::Krossover"; }
      unsigned max_production(void) { return 1; } 

      void apply(eoPopulator<t_Individual>& _pop);
      std::string print_out() const
      {
        std::ostringstream sstr;
        sstr << "KMutation rate = " << rate;
        return sstr.str();
      }
  };

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
  class Crossover : public eoGenOp<T_INDIVIDUAL> 
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIV_TRAITS t_IndivTraits;
    protected:
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      Concentration &concentration;
      Ising_CE::Structure structure;
      BitString::Crossover<t_Object> op;

    public:
      Crossover   ( Ising_CE::Structure &_str, Concentration &_c )
             : concentration(_c), structure(_str), op() {}
      Crossover   ( const Crossover &_k )
             : concentration(_k.crossover), structure(_k.structure), op(_k.op) {}
      ~Crossover() {}

      bool Load( const TiXmlElement &_node ) { return op.Load( _node ); }

      virtual std::string className() const { return op.className(); }
      unsigned max_production(void) { return 1; } 

      void apply(eoPopulator<t_Individual>& _pop)
      {
        t_Object &obj1 = (*_pop).Object();
        const t_Object &obj2 = _pop.select().Object();
        op( obj1, obj2 );
        structure << obj1;
        concentration( structure );
        obj1 << structure;
      }
      std::string print_out() const { return "TwoSites::Crossover"; }
  };
  template<class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
  class Mutation : public eoGenOp<T_INDIVIDUAL> 
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIV_TRAITS t_IndivTraits;
    protected:
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      Concentration &concentration;
      Ising_CE::Structure structure;
      BitString::Mutation<t_Object> op;

    public:
      Mutation   ( Ising_CE::Structure &_str, Concentration &_c )
             : concentration(_c), structure(_str), op() {}
      Mutation   ( const Mutation &_k )
             : concentration(_k.crossover), structure(_k.structure), op(_k.op) {}
      ~Mutation() {}

      bool Load( const TiXmlElement &_node ) { return op.Load( _node ); }

      virtual std::string className() const { return op.className(); }
      unsigned max_production(void) { return 1; } 

      void apply(eoPopulator<t_Individual>& _pop)
      {
        t_Object &obj = (*_pop).Object();
        op( obj );
        structure << obj;
        concentration( structure );
        obj << structure;
      }
      std::string print_out() const { return "TwoSites::Mutation"; }
  };

} // namespace TwoSites

#include "two_sites.impl.h"

#endif // _TWOSITES_OBJECT_H_
