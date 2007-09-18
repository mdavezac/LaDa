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

#include "gaoperators.h"

// Defines base classes for two site structures

namespace TwoSites
{
  void rearrange_structure(Ising_CE::Structure &);

  typedef SingleSite :: Object Object;

  struct Fourier
  {
    template<class T_R_IT, class T_K_IT>
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend );
    template<class T_R_IT, class T_K_IT, class T_O_IT >
    Fourier( T_R_IT _rfirst, T_R_IT _rend,
             T_K_IT _kfirst, T_K_IT _kend,
             T_O_IT _rout ); // sets rvector values from kspace values
  };

  class Concentration : public X_vs_y
  {
    public:
      types::t_real x, y;
      types::t_unsigned N;
      types::t_int Nfreeze_x, Nfreeze_y;
      std::vector<bool> sites;

    public:
      Concentration  () 
                    : X_vs_y(), x(0), y(0), N(0), Nfreeze_x(0), Nfreeze_y(0) {}
      Concentration   ( const Concentration &_conc)
                    : X_vs_y(_conc), x(_conc.x), y(_conc.y), N(_conc.N),
                      Nfreeze_x(_conc.Nfreeze_x), Nfreeze_y(_conc.Nfreeze_y),
                      sites(_conc.sites) {}
      ~Concentration() {}

      bool Load( const TiXmlElement &_node )
      {
        if( not X_vs_y::Load( _node ) ) return false;
        if( not is_singlec() )  return true;
        x = get_x();  y = get_y();
        return true;
      }

      void operator()( Ising_CE::Structure &_str );
      void operator()( Object &_obj );
      void operator()( const Ising_CE::Structure &_str, Object &_object,
                       types::t_int _concx, types::t_int _concy );
      void set( const Ising_CE::Structure &_str);

    protected:
      void normalize( Ising_CE::Structure &_str, const types::t_int _site, 
                      types::t_real _tochange);

  };

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
  class Evaluator : public GA::Evaluator< T_INDIVIDUAL, T_INDIV_TRAITS >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIV_TRAITS t_IndivTraits;
    protected:
      typedef typename t_IndivTraits::t_Object t_Object;
      typedef GA::Evaluator<t_Individual, t_IndivTraits> t_Base;
      typedef Evaluator<t_Individual, t_IndivTraits> t_This;
      typedef Traits::GAOp<t_Individual, Concentration, Fourier> t_GAOpTraits;

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
                : lattice( _c.lattice ), structure( _c.structure ),
                  lessthan(_c.lessthan), morethan(_c.moerethan),
                  concentration( _c.concentration ){}
      ~Evaluator() {};


      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Load( const TiXmlElement &_node );
      eoGenOp<t_Individual>* LoadGaOp(const TiXmlElement &_el )
       { return GA::LoadGaOp<t_GAOpTraits>( _el, structure, concentration ); }
      GA::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );

      bool initialize( t_Individual &_indiv )
      {
        GA::Random< t_GAOpTraits > random( concentration, structure, _indiv );
        _indiv.invalidate(); return true;
      }

    protected:
      bool consistency_check();
      bool Taboo(const t_Individual &_indiv );
  };



} // namespace TwoSites

#include "two_sites.impl.h"

#endif // _TWOSITES_OBJECT_H_
