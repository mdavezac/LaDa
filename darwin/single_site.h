#ifndef _SINGLE_SITE_H_
#define _SINGLE_SITE_H_

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

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif


// Defines base classes for two site structures

namespace SingleSite
{

  struct Object 
  {
    // Friend Functions
    friend void operator<<(std::string &_str, const Object &_o);
    friend void operator<<(Ising_CE::Structure &_str, const Object &_o);
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<SingleSite::Object>(SingleSite::Object &);
#endif

    // Typedefs
    public:
      typedef types::t_real t_Type;
      typedef std::vector<t_Type>  t_Container;
      typedef t_Container :: iterator iterator;
      typedef t_Container :: const_iterator const_iterator;


    // variables
    public:
      t_Container bitstring;


    // Member Functions
    public:
      Object() {}
      Object(const Object &_c) : bitstring(_c.bitstring) {};
      bool operator<<(const Ising_CE::Structure &_c)
      {
        bitstring.clear(); bitstring.reserve( _c.atoms.size() );
        Ising_CE::Structure :: t_Atoms :: const_iterator i_atom = _c.atoms.begin();
        Ising_CE::Structure :: t_Atoms :: const_iterator i_end = _c.atoms.end();
        for(; i_atom != i_end; ++i_atom )
          if ( not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T) )
            bitstring.push_back( i_atom->type > 0 ? 1.0: -1.0 );
        return true;
      }
      bool operator<<(const std::string &_c)
      {
        types::t_unsigned size = _c.size();
        bitstring.resize( size );
        std::vector<types::t_real> :: iterator i_var = bitstring.begin();
        std::vector<types::t_real> :: iterator i_end = bitstring.end();
        for(types::t_unsigned n=0; i_var != i_end; ++i_var, ++n )
          *i_var = ( _c[n] == '1' ) ? 1.0: -1.0;
        return true;
      }
      ~Object() {};
    
      bool operator==( const Object &_c ) const
      {
        return std::equal( bitstring.begin(), bitstring.end(), 
                            _c.bitstring.begin() ); 
      }
      t_Container :: const_iterator begin() const
        { return bitstring.begin();  }
      t_Container:: const_iterator end() const
        { return bitstring.end();  }
      t_Container:: iterator begin() 
        { return bitstring.begin();  }
      t_Container:: iterator end()
        { return bitstring.end();  }
      
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
      
      void print_out( std::ostream &_stream ) const
        { std::string str; str << *this; _stream << str; }
      
      void mask( types::t_unsigned _start, types::t_unsigned _end)
      {
        if ( _end > bitstring.size() )
          _end = bitstring.size();
        std::transform( bitstring.begin()+_start, bitstring.begin()+_end,
                        bitstring.begin()+_start, std::logical_not<t_Type>() );  
      }
      types::t_unsigned size()
        { return bitstring.size(); }

      t_Container& Container() { return bitstring; }
      const t_Container& Container() const { return bitstring; }
  };


  template<class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits :: Indiv<T_INDIVIDUAL> >
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
      darwin::Taboo_Base<t_Individual>* LoadTaboo(const TiXmlElement &_el );

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
  inline bool mpi::BroadCast::serialize<SingleSite::Object>( SingleSite::Object & _object )
  {
    return serialize( _object.bitstring );
  }
}
#endif

#endif // _TWOSITES_OBJECT_H_
