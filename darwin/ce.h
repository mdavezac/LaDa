#ifndef _CE_H_
#define _CE_H_

#include <string>
#include <algorithm>
#include <functional>

#include <tinyxml/tinyxml.h>
#include <eo/eoOp.h>

#include "lamarck/structure.h"
#include "lamarck/functional_builder.h"
#include "opt/opt_function_base.h"
#include "opt/types.h"

#include "evaluator.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace CE
{
  struct Object 
  {
    typedef std::vector<types::t_real> :: iterator iterator;
    typedef std::vector<types::t_real> :: const_iterator const_iterator;
    friend void operator<<(std::string &_str, const Object &_o);
    friend void operator<<(Ising_CE::Structure &_str, const Object &_o);
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<CE::Object>(CE::Object &);
#endif
    std::vector<types::t_real> bitstring;
    Object() {}
    Object(const Object &_c) : bitstring(_c.bitstring) {};
    bool operator<< (const Ising_CE::Structure &_c)
    {
      bitstring.resize( _c.atoms.size() );
      Ising_CE::Structure :: t_Atoms :: const_iterator i_atom = _c.atoms.begin();
      Ising_CE::Structure :: t_Atoms :: const_iterator i_end = _c.atoms.end();
      std::vector<types::t_real> :: iterator i_var = bitstring.begin();
      for(; i_atom != i_end; ++i_atom, ++i_var )
        *i_var = ( i_atom->type > 0 ) ? 1.0 : -1.0;
      return true;
    }
    bool operator<<(const std::string &_c)
    {
      types::t_unsigned size = _c.size();
      bitstring.resize( size );
      std::vector<types::t_real> :: iterator i_var = bitstring.begin();
      std::vector<types::t_real> :: iterator i_end = bitstring.end();
      for(types::t_unsigned n=0; i_var != i_end; ++i_var, ++n )
        *i_var = ( _c[n] == '1' ) ? 1.0 : -1.0;
      return true;
    }
    ~Object() {};
    
    types::t_real& operator[](types::t_unsigned _i) 
      { return bitstring[_i]; }
    const types::t_real& operator[](types::t_unsigned _i) const
      { return bitstring[_i]; }
    bool operator==( const Object &_c ) const
    {
      return std::equal( bitstring.begin(), bitstring.end(), 
                          _c.bitstring.begin() ); 
    }
    std::vector<types::t_real> :: const_iterator begin() const
      { return bitstring.begin();  }
    std::vector<types::t_real> :: const_iterator end() const
      { return bitstring.end();  }
    std::vector<types::t_real> :: iterator begin() 
      { return bitstring.begin();  }
    std::vector<types::t_real> :: iterator end()
      { return bitstring.end();  }

    types::t_real get_concentration() const
    {
      std::vector<types::t_real> :: const_iterator i_var = bitstring.begin();
      std::vector<types::t_real> :: const_iterator i_var_end = bitstring.end();
      types::t_real result = 0.0;
      for(; i_var != i_var_end; ++i_var )
        result += *i_var;
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
                      bitstring.begin()+_start, std::negate<types::t_real>() );  
    }
    types::t_unsigned size()
      { return bitstring.size(); }
  };


  class Evaluator : public darwin::Evaluator<Object>, public VA_CE :: Functional_Builder
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<CE::Evaluator>( CE::Evaluator & );
#endif
    public:
      typedef function::Base<> t_Functional;
      // typedef t_Object is inherited from darwin::Evaluator<...> 
      
    protected:
      using VA_CE::Functional_Builder::t_VA_Functional;
    public:
      using darwin::Evaluator<Object> :: Load;

    protected:
      t_VA_Functional functional;
      types::t_real crossover_probability;
      Ising_CE::Structure structure;
      bool single_concentration;
      types::t_real x;
      types::t_real lessthan, morethan;

    public:
      Evaluator() : functional(), crossover_probability(0.5),
                    single_concentration(false), x(0), lessthan(-1.0), morethan(1.0) {}; 
      ~Evaluator() 
      {
        if ( functional.get_functional1() ) 
          delete functional.get_functional1();
        if ( functional.get_functional2() ) 
          delete functional.get_functional2();
      };

      void* const init( t_Object &_object );
      bool Load( const TiXmlElement &_node );
      bool Load ( t_Object &_indiv, const TiXmlElement &_node, bool _type );
      bool Save ( const t_Object &_indiv, TiXmlElement &_node, bool _type ) const;
      void LoadAttribute ( const TiXmlAttribute &_att );
      eoOp<t_Object>* LoadGaOp(const TiXmlElement &_el );
      eoMonOp<const Object>* LoadTaboo(const TiXmlElement &_el );

      bool Krossover( t_Object  &_offspring, const t_Object &_parent,
                      bool _range = false );
      bool Crossover( t_Object &_obj1, const t_Object &_obj2 );

      bool initialize( t_Object &_object );
      void set_object( t_Object &_object, const void * const _f ) {}
      void* const LoadMinimizer(const TiXmlElement &_el );
      bool Taboo(const t_Object &_object );
  };

} // namespace CE


#endif // _CE_OBJECT_H_
