#include <functional>
#include <algorithm>
#include <ext/algorithm>

#include "ce.h"
#include "functors.h"
#include "print_xmgrace.h"
#include <lamarck/atom.h>
#include <opt/va_minimizer.h>

namespace CE
{
  void operator<<(std::string &_str, const Object &_o)
  {
    std::vector<types::t_real> :: const_iterator i_var = _o.bitstring.begin();
    std::vector<types::t_real> :: const_iterator i_end = _o.bitstring.end();
    std::ostringstream sstr;
    for(; i_var != i_end; ++i_var )
      sstr << ( ( *i_var >  0.0 ) ? '1' : '0' );
    _str = sstr.str();
  }
  void operator<<(Ising_CE::Structure &_str, const Object &_o)
  {
    std::vector<types::t_real> :: const_iterator i_var = _o.bitstring.begin();
    std::vector<types::t_real> :: const_iterator i_end = _o.bitstring.end();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure :: t_Atoms :: iterator i_atom_end = _str.atoms.end();
    for(; i_var != i_end and i_atom != i_atom_end; ++i_var, ++i_atom )
      i_atom->type = ( *i_var > 0 ) ? 1.0 : -1.0;
  }

  bool Evaluator :: Load( Object &_obj, const TiXmlElement &_node, bool _type )
  {
    if ( _type == darwin::LOADSAVE_SHORT )
    {
      if( not _node.Attribute("string") )
        return false;
      _obj = std::string(_node.Attribute("string"));
      return true;
    }

    Ising_CE::Structure s; 
    if ( not s.Load(_node) )
      return false;
    _obj = s;
    return true;
  }
  bool Evaluator :: Save( const Object &_obj, TiXmlElement &_node, bool _type ) const
  {
    if ( _type == darwin::LOADSAVE_SHORT )
    {
      std::string str; str << _obj;
      _node.SetAttribute("string", str.c_str());
      return true;
    }

    Ising_CE::Structure s = structure; 
    s << _obj;
    Ising_CE::fourrier_to_kspace( s.atoms.begin(),  s.atoms.end(),
                                  s.k_vecs.begin(), s.k_vecs.end() );
    s.print_xml(_node);

    return true;
  }
  bool Evaluator :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *functional_xml = _node.FirstChildElement("Functional");
    for(; functional_xml; functional_xml = functional_xml->NextSiblingElement("Functional") )
    {
      std::string str = ""; 
      if ( functional_xml->Attribute("type") )
        str = functional_xml->Attribute("type");
      if ( str.compare("CE") == 0 )
        break;
    }
    if ( not functional_xml )
    {
      std::cerr << "No <Functional type=\"CE\"> found "
                << std::endl << "Giving up" << std::endl;
      return false;
    }

    if ( not VA_CE::Functional_Builder::Load(*functional_xml) )
      return false;

    // reads structure from input
    const TiXmlElement *structure_xml = _node.FirstChildElement( "Structure" );
    if ( not structure_xml )
      structure_xml = functional_xml->FirstChildElement( "Structure" );
    if ( not structure_xml )
    {
      std::cerr << "Could not find Structure " << std::endl;
      return false;
    }
    structure.Load(*structure_xml);

    add_equivalent_clusters();

    generate_functional(structure, &functional);
    structure.set_kvectors( functional.get_functional2()->get_kvectors() );

    return true;
  }

  bool Evaluator :: Crossover ( t_Object &_obj1, const t_Object &_obj2 )
  {
    std::vector<types::t_real> :: iterator i_var1 = _obj1.begin();
    std::vector<types::t_real> :: const_iterator i_var2 = _obj2.begin();
    std::vector<types::t_real> :: const_iterator i_var2_end = _obj2.end();
    for(; i_var2 != i_var2_end; ++i_var1, ++i_var2)
      if ( rng.flip(crossover_probability) ) 
        *i_var1 = *i_var2;
    return true;
  }

  // expects kspace value to exist!!
  bool Evaluator :: Krossover( Object  &_offspring, const Object &_parent,
                               bool _range )
  {
    Ising_CE::Structure str1 = structure, str2 = structure;
    str1 << _offspring; str2 << _parent;
    Ising_CE::fourrier_to_kspace( str1.atoms.begin(),  str1.atoms.end(),
                                  str1.k_vecs.begin(), str1.k_vecs.end() );
    Ising_CE::fourrier_to_kspace( str2.atoms.begin(),  str2.atoms.end(),
                                  str2.k_vecs.begin(), str2.k_vecs.end() );
    do 
    {
      if ( _range and str1.k_vecs.size() > 2 ) // range crossover ... kvec should be oredered according to size
      {  
        types::t_unsigned n = (types::t_unsigned)
           std::floor(   (types::t_real) rng.random ( str1.k_vecs.size() - 1 ) 
                       * (types::t_real) crossover_probability );
        __gnu_cxx::copy_n( str2.k_vecs.begin(), n, str1.k_vecs.begin() );
      }
      else // every point crossover
      {
        Ising_CE::Structure::t_kAtoms :: const_iterator i_p = str2.k_vecs.begin();
        Ising_CE::Structure::t_kAtoms :: const_iterator i_p_end = str2.k_vecs.end();
        Ising_CE::Structure::t_kAtoms :: iterator i_o = str1.k_vecs.begin();
        for ( ; i_p != i_p_end; ++i_p, ++i_o)
          if ( rng.flip(crossover_probability) ) 
            i_o->type = i_p->type;
      }
  
      if ( set_concentration( str1 ) )
        break;
      str1 << _offspring;
    }
    while (true);
    _offspring = str1;

    return true; // offspring has changed!
  }

  bool Evaluator :: set_concentration( Ising_CE::Structure &_object )
  {
    typedef std::complex<types::t_real> t_complex;
    types::t_real N = 1.0 / types::t_real( _object.atoms.size() );
    types::t_real tolerance = N;
    types::t_unsigned count = 0;

    // Then FT back to r space
    t_complex  *Chold = new t_complex[ _object.atoms.size() ];
    t_complex  *i_chold = Chold;
    Ising_CE::fourrier_to_rspace( _object.atoms.begin(), _object.atoms.end(),
                                  _object.k_vecs.begin(), _object.k_vecs.end(),
                                  i_chold );

    // Finally, normalizes and copies
    i_chold = Chold;
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _object.atoms.end();
    Ising_CE::Structure::t_Atoms::iterator i_atom_begin = _object.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom = i_atom_begin;
    types::t_real X = 0;
    for (; i_atom != i_atom_end; ++i_atom, ++i_chold)
    {
      i_atom->type = ( std::real(*i_chold) > 0 ) ? 1.0 : -1.0;
      X += i_atom->type;
    }
    X *= N;
    if (    ( not single_concentration )
         or ( std::abs(X - x) < tolerance ) )
      return true;

    bool sign = (X > x ); 
    types::t_real step = ( sign ? -1: 1 );
    while( ( sign and X > x ) or ( (not sign) and X < x ) )
    {
      X = 0; i_atom = i_atom_begin;
      _object.k_vecs.front().type += step;
      i_chold = Chold; i_atom = i_atom_begin; X = 0.0;
      Ising_CE::fourrier_to_rspace( _object.atoms.begin(), _object.atoms.end(),
                                    _object.k_vecs.begin(), _object.k_vecs.end(),
                                    i_chold );
      for (; i_atom != i_atom_end; ++i_atom, ++i_chold)
        X += ( std::real(*i_chold) > 0 ) ? 1.0 : -1.0;
      X *= N;
    }
    while ( std::abs(X - x) > tolerance and count < 64 )
    {
      ++ count;
      step = ( X > x ) ? -0.5 * std::abs(step) : 0.5 * std::abs(step);
      _object.k_vecs.front().type += step;
      i_chold = Chold; i_atom = i_atom_begin; X = 0.0;
      Ising_CE::fourrier_to_rspace( _object.atoms.begin(), _object.atoms.end(),
                                    _object.k_vecs.begin(), _object.k_vecs.end(),
                                    i_chold );
      for (; i_atom != i_atom_end; ++i_atom, ++i_chold)
        X += ( std::real(*i_chold) > 0 ) ? 1.0 : -1.0;
      X *= N;
    }

    i_chold = Chold; i_atom = i_atom_begin; X = 0;
    for (; i_atom != i_atom_end; ++i_atom, ++i_chold)
    {
      i_atom->type = ( std::real(*i_chold) > 0 ) ? 1.0 : -1.0;
      X += ( std::real(*i_chold) > 0 ) ? 1.0 : -1.0;
    }
    X *= N;
    delete[] Chold;
    return std::abs(X - x) < tolerance;
  }

  eoOp<Object>* Evaluator::LoadGaOp(const TiXmlElement &_el )
  {
    std::string value = _el.Value();

    if ( value.compare("Crossover") == 0 )
    {
      _el.Attribute( "prob", &crossover_probability );
      crossover_probability = crossover_probability > 0 ? std::abs(crossover_probability) : 0.5;
      if ( crossover_probability > 1 ) crossover_probability = 0.5;
      std::ostringstream sstr;
      sstr << "Crossover rate = " << crossover_probability;
      darwin::printxmg.add_comment(sstr.str());
      // pointer is owned by caller !!
      return new darwin::mem_binop_t<Evaluator, Object, void>( *this, &Evaluator::Crossover, 
                                                               std::string( "Crossover" ) );
    }
    else if ( value.compare( "Krossover" ) == 0 )
    {
      bool att = false;
      _el.Attribute( "prob", &crossover_probability );
      crossover_probability = crossover_probability > 0 ? std::abs(crossover_probability) : 0.5;
      if ( crossover_probability > 1 ) crossover_probability = 0.5;
      std::ostringstream sstr;
      sstr << "Krossover, rate = " << crossover_probability;
      darwin::printxmg.add_comment(sstr.str());
      if ( _el.Attribute("type") )
      {
        std::string str =  _el.Attribute("type");
        if ( str.compare("range") == 0 ) 
          { att = true; darwin::printxmg.add_to_last( ", Range = true" ); }
      }
      // pointer is owned by caller !!
      return new darwin::mem_binop_t<Evaluator, Object, bool>( *this, &Evaluator::Krossover, 
                                                               std::string( "Krossover" ), att);
    }

    return NULL;
  }
  void Evaluator :: LoadAttribute ( const TiXmlAttribute &_att )
  {
    std::string str = _att.Name();
    if ( str.compare("x") == 0 )
    {
      x = 2.0 * _att.DoubleValue() - 1.0;
      single_concentration = true;
      if(    std::abs(x-1.0) < types::tolerance 
          or std::abs(x-1.0) < types::tolerance  ) 
      {
        std::cerr << " wrong concentration  on input x = "
                  << _att.DoubleValue() << std::endl;
        exit(0);
      }
    }
  }


  bool Evaluator::initialize( Object &_object )
  {
    Ising_CE::Structure str = structure;
    _object.bitstring.resize( structure.atoms.size() );
    do
    {
     Object :: iterator i_var_end = _object.end();
     Object :: iterator i_var = _object.begin();
     for(; i_var != i_var_end; ++i_var )
       *i_var = rng.flip() ? 1.0 : -1.0;
     if ( not single_concentration )
       return true;
     str << _object;
     Ising_CE::fourrier_to_kspace( str.atoms.begin(),  str.atoms.end(), 
                                   str.k_vecs.begin(), str.k_vecs.end() );
    } while(not set_concentration( str ) );
    return true;
  }

  void* const Evaluator::init( Object &_object )
  {
    functional.set_variables( &_object.bitstring );
    return &functional;
  }

  void* const Evaluator::LoadMinimizer(const TiXmlElement &_el )
  {
    std::string name = _el.Value();
    if ( name.compare("Minimizer") )
      return NULL;

    if ( not _el.Attribute("type") )
      return false;
    
    name = _el.Attribute("type");

    if ( name.compare("VA") == 0 )
    {
      darwin::printxmg.add_comment("VA optimizer");
      // pointer is owned by caller !!
      // don't deallocate
      return new minimizer::VA<t_Functional>( _el );
    }
    else if ( name.compare("SA") == 0 )
    {
      darwin::printxmg.add_comment("SA optimizer");
      // pointer is owned by caller !!
      // don't deallocate
      return new minimizer::VA<t_Functional>( _el );
    }
    else if ( name.compare("Beratan") == 0 )
    {
      darwin::printxmg.add_comment("SA optimizer");
      // pointer is owned by caller !!
      // don't deallocate
      return new minimizer::Beratan<t_Functional>( _el );
    }
    

    return NULL;
  }

} // namespace CE



#ifdef _MPI
namespace mpi
{
  template<>
  bool mpi::BroadCast::serialize<CE::Evaluator>( CE::Evaluator & _ev )
  {
    if( not serialize( _ev.crossover_probability ) ) return false;
    if( not serialize( _ev.structure ) ) return false;
    if( not serialize( _ev.single_concentration ) ) return false;
    if( not serialize( _ev.x ) ) return false;
    if( not serialize<VA_CE::Functional_Builder>(_ev) ) return false;

    if ( stage == COPYING_FROM_HERE and rank() != ROOT_NODE )
    {
      _ev.generate_functional( _ev.structure, &_ev.functional );
      _ev.structure.set_kvectors( _ev.functional.get_functional2()->get_kvectors() );
    }

    return true;
  }
  template<>
  bool mpi::BroadCast::serialize<CE::Object>( CE::Object & _object )
  {
    return serialize( _object.bitstring );
  }
}
#endif
