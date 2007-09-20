//
//  Version: $Id$
//

#ifndef _GAOPERATORS_H_
#define _GAOPERATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <algorithm>
#include <ext/algorithm>
#include <functional>
#include <string>
#include <sstream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include "lamarck/structure.h"
#include "opt/types.h"

#include "functors.h"
#include "gatraits.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace GA
{
  template<class T_GAOPTRAITS> 
  class Krossover : public eoGenOp<typename T_GAOPTRAITS :: t_Individual> 
  {
    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_GAOpTraits :: t_FourierRtoK t_FourierRtoK;
      typedef typename t_GAOpTraits :: t_FourierKtoR t_FourierKtoR;
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      t_Concentration &concentration;
      Ising_CE::Structure &structure;
      types::t_real rate;
      bool do_range;

    public:
      Krossover   ( t_Concentration &_c, Ising_CE::Structure &_str )
                : concentration(_c), structure(_str), rate(0.5), do_range(false) {}
      Krossover   ( const Krossover &_k )
                : concentration(_k.concentration), structure(_k.structure),
                  rate(_k.rate), do_range(_k.do_range) {}
      ~Krossover() {}

      bool Load( const TiXmlElement &_node );

      virtual std::string className() const { return "TwoSites::Krossover"; }
      unsigned max_production(void) { return 1; } 

      bool operator()( t_Individual &_indiv, const t_Individual &_parent );
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop, _pop.select() ) ) (*_pop).invalidate(); }
      std::string print_out() const
      {
        std::ostringstream sstr;
        sstr << "Krossover, rate = " << rate;
        if (do_range) sstr << ", Range = true ";
        return sstr.str();
      }
  };
  template<class T_GAOPTRAITS>
  bool Krossover<T_GAOPTRAITS> :: Load( const TiXmlElement &_node )
  {
    std::string name = _node.Value();
    if ( name.compare("Krossover") ) return false;

    _node.Attribute( "prob", &rate );
    rate = rate > 0 ? std::abs(rate) : 0.5;
    if ( rate > 1 ) rate = 0.5;
    if ( _node.Attribute("type") )
    {
      std::string str =  _node.Attribute("type");
      if ( str.compare("range") == 0 )  do_range = true;
    }
    return true;
  }
  // expects kspace value to exist!!
  template<class T_GAOPTRAITS>
  bool Krossover<T_GAOPTRAITS> :: operator()( t_Individual &_indiv, const t_Individual &_parent )
  {
    t_Object &offspring  = _indiv.Object();
    const t_Object &parent  = _parent.Object();
    Ising_CE::Structure str1 = structure, str2 = structure;
    str1 << offspring; str2 << parent;
    t_FourierRtoK( str1.atoms.begin(),  str1.atoms.end(),
                   str1.k_vecs.begin(), str1.k_vecs.end() );
    t_FourierKtoR( str2.atoms.begin(),  str2.atoms.end(),
                   str2.k_vecs.begin(), str2.k_vecs.end() );

    // range crossover ... kvec should be oredered according to size
    if ( do_range and str1.k_vecs.size() > 2 ) 
    {  
      types::t_unsigned n = (types::t_unsigned)
         std::floor(   (types::t_real) rng.random ( str1.k_vecs.size() - 1 ) 
                     * (types::t_real) rate );
      __gnu_cxx::copy_n( str2.k_vecs.begin(), n, str1.k_vecs.begin() );
    }
    else // every point crossover
    {
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p = str2.k_vecs.begin();
      Ising_CE::Structure::t_kAtoms :: const_iterator i_p_end = str2.k_vecs.end();
      Ising_CE::Structure::t_kAtoms :: iterator i_o = str1.k_vecs.begin();
      for ( ; i_p != i_p_end; ++i_p, ++i_o)
        if ( rng.flip(rate) ) 
          i_o->type = i_p->type;
    }
  
    concentration( str1 );
    offspring << str1;

    return true;
  }

  template<class T_GAOPTRAITS>
  class KMutation : public eoGenOp<typename T_GAOPTRAITS::t_Individual> 
  {
    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_GAOpTraits :: t_FourierRtoK t_FourierRtoK;
      typedef typename t_GAOpTraits :: t_FourierKtoR t_FourierKtoR;
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      t_Concentration &concentration;
      Ising_CE::Structure &structure;
      types::t_real rate;

    public:
      KMutation   ( t_Concentration &_c, Ising_CE::Structure &_str )
                : concentration(_c), structure(_str), rate(0.5) {}
      KMutation   ( const KMutation &_k )
                : concentration(_k.concentration), structure(_k.structure),
                  rate(_k.rate) {}
      ~KMutation() {}

      bool Load( const TiXmlElement &_node );

      virtual std::string className() const { return "TwoSites::KMutation"; }
      unsigned max_production(void) { return 1; } 

      bool operator()( t_Individual &_indiv );
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }
      std::string print_out() const
      {
        std::ostringstream sstr;
        sstr << "KMutation, rate = " << rate;
        return sstr.str();
      }
  };
  template<class T_GAOPTRAITS> 
  bool KMutation<T_GAOPTRAITS> :: Load( const TiXmlElement &_node )
  {
    std::string name = _node.Value();
    if ( name.compare("KMutation") ) return false;

    _node.Attribute( "prob", &rate );
    rate = rate > 0 ? std::abs(rate) : 0.5;
    if ( rate > 1 ) rate = 0.5;
    return true;
  }
  // expects kspace value to exist!!
  template<class T_GAOPTRAITS>
  bool KMutation<T_GAOPTRAITS> :: operator()(t_Individual &_indiv )
  {
    t_Object &object  = _indiv.Object();
    Ising_CE::Structure str = structure;
    str << object;
    t_FourierRtoK( str.atoms.begin(),  str.atoms.end(),
                   str.k_vecs.begin(), str.k_vecs.end() );

    
    Ising_CE::Structure::t_kAtoms :: iterator i_k = str.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms :: iterator i_k_end = str.k_vecs.end();
    types::t_real max = std::norm( i_k->type );
    for(; i_k != i_k_end; ++i_k)
      if ( max < std::norm( i_k->type ) ) max = std::norm( i_k->type );

    i_k = str.k_vecs.begin();
    max = std::sqrt( max );
    bool result = false;
    for(; i_k != i_k_end; ++i_k)
      if ( rng.flip(rate) )
      {
        i_k->type = std::complex<types::t_real>( rng.uniform( 2.0 * max ) - max, 
                                                 rng.uniform( 2.0 * max ) - max );
        result = true;
      }
  
    concentration( str );
    object << str;
    return true;
  }
  template<class T_GAOPTRAITS> 
  class KRandom : public eoGenOp<typename T_GAOPTRAITS :: t_Individual> 
  {
    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      t_Concentration &concentration;
      Ising_CE::Structure structure;

    public:
      KRandom   ( t_Concentration &_c, Ising_CE::Structure &_str )
              : concentration(_c), structure(_str) {}
      KRandom   ( t_Concentration &_c, Ising_CE::Structure &_str, t_Individual &_indiv )
              : concentration(_c), structure(_str) { operator()(_indiv); }
      KRandom   ( const KRandom &_k )
              : concentration(_k.crossover), structure(_k.structure) {}
      ~KRandom() {}

      bool Load( const TiXmlElement &_node ) { return true; }

      virtual std::string className() const { return "Darwin::Random"; }
      unsigned max_production(void) { return 1; } 

      bool operator()( t_Individual &_indiv )
      {
        t_Object &obj = _indiv.Object();

        Ising_CE::Structure::t_Atoms :: const_iterator i_kvec = structure.k_vecs.begin();
        Ising_CE::Structure::t_Atoms :: const_iterator i_kvec_end = structure.k_vecs.end();
        types::t_real n = 3.0 / (types::t_real) structure.k_vecs.size(); 
        for(; i_kvec != i_kvec_end; ++i_kvec )
          if( rng.flip(n) )
            i_kvec->type = std::complex<types::t_real>( rng.uniform( 5.0 ) - 2.5,
                                                        rng.uniform( 5.0 ) - 2.5  );

        concentration( structure );
        obj << structure;

        return true;
      }

      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }
      std::string print_out() const { return "Darwin::KRandom"; }
  };



  template<class T_GAOPTRAITS> 
  class Crossover : public eoGenOp<typename T_GAOPTRAITS :: t_Individual > 
  {

    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_Object t_Object;
    protected:
      t_Concentration &concentration;
      Ising_CE::Structure structure;
      BitString::Crossover<t_Object> op;

    public:
      Crossover   ( t_Concentration &_c, Ising_CE::Structure &_str )
             : concentration(_c), structure(_str), op() {}
      Crossover   ( const Crossover &_k )
             : concentration(_k.crossover), structure(_k.structure), op(_k.op) {}
      ~Crossover() {}

      bool Load( const TiXmlElement &_node ) { return op.Load( _node ); }

      virtual std::string className() const { return op.className(); }
      unsigned max_production(void) { return 1; } 

      bool operator()(t_Individual& _indiv, const t_Individual _parent )
      {
        t_Object &obj1 = _indiv.Object();
        const t_Object &obj2 = _parent.Object();
        op( obj1, obj2 );
        concentration( obj1 );
        return true;
      }
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop, _pop.select() ) ) (*_pop).invalidate(); }
      std::string print_out() const
      {
        std::ostringstream sstr;
        sstr << "Crossover, rate = " << op.rate;
        return sstr.str();
      }
  };
  template<class T_GAOPTRAITS> 
  class Mutation : public eoGenOp<typename T_GAOPTRAITS :: t_Individual> 
  {
    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      t_Concentration &concentration;
      Ising_CE::Structure structure;
      BitString::Mutation<t_Object> op;

    public:
      Mutation   ( t_Concentration &_c, Ising_CE::Structure &_str )
             : concentration(_c), structure(_str), op() {}
      Mutation   ( const Mutation &_k )
             : concentration(_k.crossover), structure(_k.structure), op(_k.op) {}
      ~Mutation() {}

      bool Load( const TiXmlElement &_node ) { return op.Load( _node ); }

      virtual std::string className() const { return op.className(); }
      unsigned max_production(void) { return 1; } 

      bool operator()( t_Individual &_indiv )
      {
        t_Object &obj = _indiv.Object();
        op( obj );
        concentration( obj );
        return true;
      }
      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }
      std::string print_out() const 
      {
        std::ostringstream sstr;
        sstr << "Mutation, rate = " << op.rate;
        return sstr.str();
      }
  };
  template<class T_GAOPTRAITS> 
  class Random : public eoGenOp<typename T_GAOPTRAITS :: t_Individual> 
  {
    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_Object t_Object;

    protected:
      t_Concentration &concentration;
      Ising_CE::Structure structure;

    public:
      Random   ( t_Concentration &_c, Ising_CE::Structure &_str )
             : concentration(_c), structure(_str) {}
      Random   ( t_Concentration &_c, Ising_CE::Structure &_str, t_Individual &_indiv )
             : concentration(_c), structure(_str) { operator()(_indiv); }
      Random   ( const Random &_k )
             : concentration(_k.crossover), structure(_k.structure) {}
      ~Random() {}

      bool Load( const TiXmlElement &_node ) { return true; }

      virtual std::string className() const { return "Darwin::Random"; }
      unsigned max_production(void) { return 1; } 

      bool operator()( t_Individual &_indiv )
      {
        t_Object &obj = _indiv.Object();
        obj.bitstring.clear();

        Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
        Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom )
          if ( not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T) )
            obj.bitstring.push_back( rng.flip() ? 1.0: -1.0 );

        concentration( obj );

        return true;
      }

      void apply(eoPopulator<t_Individual>& _pop)
        { if ( operator()( *_pop ) ) (*_pop).invalidate(); }
      std::string print_out() const { return "Darwin::Random"; }
  };

  template<class T_GAOPTRAITS>
    eoGenOp<typename T_GAOPTRAITS::t_Individual >*
      LoadGaOp(const TiXmlElement &_el, Ising_CE::Structure &_structure, 
               typename T_GAOPTRAITS :: t_Concentration &_concentration )
  {
    std::string value = _el.Value();

    if ( value.compare("Crossover") == 0 )
    {
      Crossover<T_GAOPTRAITS> *crossover
         = new Crossover<T_GAOPTRAITS>(_concentration, _structure);
      if ( not crossover ) goto errorout;
      if ( not crossover->Load( _el ) ) 
      {
        delete crossover;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << crossover->print_out() << Print::endl;
      // pointer is owned by caller !!
      return crossover;
    }
    else if ( value.compare("Mutation") == 0 )
    {
      Mutation<T_GAOPTRAITS> *mutation
         = new Mutation<T_GAOPTRAITS>(_concentration, _structure);
      if ( not mutation ) goto errorout;
      if ( not mutation->Load( _el ) ) 
      {
        delete mutation;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << mutation->print_out() << Print::endl;
      // pointer is owned by caller !!
      return mutation;
    }
    else if ( value.compare( "Random" ) == 0 )
    {
      Random<T_GAOPTRAITS> *random
         = new Random<T_GAOPTRAITS>(_concentration, _structure);
      Print::xmg << Print::Xmg::comment << random->print_out() << Print::endl;
      // pointer is owned by caller !!
      return random;
    }
    else if ( value.compare( "Krossover" ) == 0 )
    {
      Krossover<T_GAOPTRAITS> *krossover
         = new Krossover<T_GAOPTRAITS>( _concentration, _structure );
      if ( not krossover ) goto errorout;
      if ( not krossover->Load( _el ) ) 
      {
        delete krossover;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << krossover->print_out() << Print::endl;
      // pointer is owned by caller !!
      return krossover;
    }
    else if ( value.compare( "KMutation" ) == 0 )
    {
      KMutation<T_GAOPTRAITS> *kmutation
         = new KMutation<T_GAOPTRAITS>( _concentration, _structure );
      if ( not kmutation ) goto errorout;
      if ( not kmutation->Load( _el ) ) 
      {
        delete kmutation;
        return NULL;
      }
      Print::xmg << Print::Xmg::comment << kmutation->print_out() << Print::endl;
      // pointer is owned by caller !!
      return kmutation;
    }

    return NULL;

errorout:
    std::cerr << " Memory Allocation error while creating GA Operators " << std::endl;
    return NULL;
  }



  template< class T_GAOPTRAITS >
  class xTaboo : public Taboo_Base< typename T_GAOPTRAITS::t_Individual >
  {
    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_GAOpTraits :: t_Object t_Object;

    protected:
      t_Concentration &concentration;
      types::t_real morethan, lessthan;

    public:
       xTaboo   ( t_Concentration &_c) 
              : concentration(_c) {}
      ~xTaboo() {}

      std::string className() const { return "Darwin::xTaboo"; }

      bool operator()( const t_Individual& _indiv ) const
      {
        concentration.set( _indiv.Object() );
        return concentration.x < lessthan and concentration.x > morethan;
      }
      bool Load( const TiXmlElement &_el ) 
      {
        const TiXmlElement *child = &_el;
        std::string name = _el.Value();
        if ( name.compare("Concentration") )
          child = _el.FirstChildElement( "Concentration" );
        if ( not child ) return false;
        
        double d; 
        if ( child->Attribute( "lessthan" ) )
          child->Attribute( "xlessthan", &d );
        lessthan = ( d > 0 and d < 1 ) ? 2.0*d-1.0: 1.0;
        if ( child->Attribute( "morethan" ) )
          child->Attribute( "morethan", &d );
        morethan = ( d > 0 and d < 1 ) ? 2.0*d-1.0: -1.0;
        if ( lessthan < morethan ) return false;
       
        Print::xmg << Print::Xmg::comment << Print::fixed << Print::setprecision(3) 
                   << "Taboo x in [ " << 0.5*(morethan+1.0)
                   << ", "  << 0.5*(lessthan+1.0) << "] " << Print::endl;
        return true;
      }
  };
}
#endif
