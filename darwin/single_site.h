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
#include "concentration.h"
#include "functors.h"
#include "taboos.h"
#include "bitstring.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif


// Defines base classes for two site structures

namespace SingleSite
{

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
    
      
      types::t_real get_concentration() const
      {
        t_Container :: const_iterator i_var = bitstring.begin();
        t_Container :: const_iterator i_var_end = bitstring.end();
        types::t_real result = 0.0;
        for(; i_var != i_var_end; ++i_var )
          result += *i_var > 0 ? 1.0: -1.0;
        result /= static_cast<types::t_real>(bitstring.size());
        return result;
      }
  };

  std::ostream& operator<<(std::ostream &_stream, const Object &_o);
  void operator<<(std::string &_str, const Object &_o);
  void operator<<(Object &_o, const std::string &_str );
  void operator<<(Ising_CE::Structure &_str, const Object &_o);
  void operator<<(Object &_o, const Ising_CE::Structure &_c);
  void operator<<(Object &_o, const Ising_CE::Structure &_str);

  template<class T_INDIVIDUAL>
  class Evaluator : public GA::Evaluator< T_INDIVIDUAL >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
    protected:
      typedef typename t_Individual::t_IndivTraits t_IndivTraits;
      typedef typename t_IndivTraits::t_Object t_Object;
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
      types::t_real crossover_probability;
      types::t_real x;
      types::t_real lessthan, morethan;
      bool singlec;

    public:
      Evaluator   ()
                : crossover_probability(0.5), 
                  x(0), lessthan(1.0), morethan(-1.0), singlec(false) {}
      Evaluator   ( const t_Base &_c )
                : crossover_probability( &_c.crossover_probability ), 
                  x(_c.x), lessthan(_c.lessthan), morethan(_c.moerethan), singlec(_c.singlec) {}
      ~Evaluator() {};


      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el );
      GA::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );

      bool Krossover( t_Individual  &_offspring, const t_Individual &_parent,
                      bool _range = false );
      bool Crossover( t_Individual &_indiv1, const t_Individual &_indiv2 );

      bool initialize( t_Individual &_indiv );

    protected:
      bool consistency_check();
      void set_concentration( Ising_CE::Structure &_str );
      void normalize( Ising_CE::Structure &_str, 
                      types::t_real _tochange);
      bool Taboo(const t_Individual &_indiv );
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
