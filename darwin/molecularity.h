//
//  Version: $Id$
//
#ifndef _MOLECULARITY_H_
#define _MOLECULARITY_H_

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

#include "pescan.h"
#include "evaluator.h"
#include "concentration.h"
#include "functors.h"
#include "individual.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace Molecularity
{
  // Exact same object as BandGap, except for the two friends
  class Object : public BitString::Object<>, Vff::Keeper, Pescan::Keeper
  {
    protected:
      typedef BitString :: Object<> t_Base;

    public:
      typedef t_Base :: t_Type t_Type;
      typedef t_Base :: t_Container t_Container;

    public:
      Object() {}
      Object(const Object &_c) : t_Base(_c), Vff::Keeper(_c), Pescan::Keeper(_c) {};
      ~Object() {};
      
      types::t_int get_concentration() const
      {
        t_Container :: const_iterator i_var = bitstring.begin();
        t_Container :: const_iterator i_var_end = bitstring.end();
        types::t_int result = 0;
        for(; i_var != i_var_end; ++i_var )
          result += *i_var > 0 ? 1: -1;
        return result;
      }

      bool Load( const TiXmlElement &_node )
        { return Vff::Keeper::Load(_node) and Pescan::Keeper::Load(_node); }
      bool Save( TiXmlElement &_node ) const
        { return Vff::Keeper::Save(_node) and Pescan::Keeper::Save(_node); }
  };

  class Concentration
  {
    public:
      types::t_real x0;
      bool single_c;
      types::t_real x;

    public:
      Concentration () : x0(-1), single_c(false) {}

      types::t_real operator()( Ising_CE::Structure );
      void set( Ising_CE::Structure );
      bool Load ( const TiXmlElement &_node );

    protected:
      void normalize( Ising_CE::Structure &_str, types::t_real _tochange );
  }

  std::ostream& operator<<(std::ostream &_stream, const Object &_o);
  void operator<<(std::string &_str, const Object &_o);
  void operator<<(Object &_o, const std::string &_str );
  void operator<<(Ising_CE::Structure &_str, const Object &_o);
  void operator<<(Object &_o, const Ising_CE::Structure &_c);
  void operator<<(Object &_o, const Ising_CE::Structure &_str);


  class Evaluator : public Darwin::Evaluator< Individual::Types<Object>::Vector >
  {
    protected:
      typedef Darwin::Evaluator< Individual::Types<Object>::Vector > t_Base;
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;
    public:
      typedef  t_Base :: t_Individual t_Individual;

    public:
      using t_Base :: Load;
      using t_Base :: Save;
    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      Ising_CE::Lattice lattice;
      Ising_CE::Structure structure;
      Concentration concentration;
      Pescan::Darwin pescan;
      Pescan::Vff vff;

    public:
      Evaluator() : t_Base(), pescan(_structure), vff(_structure) {}
      Evaluator   ( const t_Base &_c )
                : lattice( _c.lattice ), structure( _c.structure ),
                  lessthan(_c.lessthan), morethan(_c.moerethan), 
                  concentration( _c.concentration ), vff(_c.vff), pescan(_c.pescan) {}
      ~Evaluator() {}

      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el );
      darwin::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );

      bool initialize( t_Individual &_indiv );

    protected:
      bool consistency_check();
  }

} // namespace BandGap


#include "molecularity.impl.h"

#ifdef _MPI
namespace mpi
{
  template<>
  inline bool mpi::BroadCast::serialize<BandGap::Object>( BandGap::Object & _object )
  {
    return     serialize<Pescan::Keeper>( _object )
           and serialize<Pescan::Vff>( _object )
           and _object.broadcast( *this );
  }
}
#endif

#endif // _MOLECULARITY_H_
