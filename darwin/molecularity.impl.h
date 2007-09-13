
namespace Molecularity
{

  template<class T_INDIVIDUAL, class T_INDIV_TRAITS>
  bool Evaluator<T_INDIVIDUAL,T_INDIV_TRAITS>::initialize( t_Individual &_indiv )
  {
    t_Object &object = _indiv.Object();
    object.bitstring.clear(); 
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom = structure.atoms.begin();
    Ising_CE::Structure::t_Atoms :: const_iterator i_atom_end = structure.atoms.end();
    types::t_int concx = 0;
    for(; i_atom != i_atom_end; ++i_atom )
    {
      bool flip = rng.flip();
      if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
        flip = ( i_atom->type > 0 );
      object.bitstring.push_back( flip ? -1.0: 1.0 );
      flip ? ++concx: --concx;
    }
    if ( singlec )
    {
      types::t_unsigned N = structure.atoms.size() >> 1; 
      types::t_real xto_change = (types::t_real) N * x  - concx;
      if ( xto_change > -1.0 and xto_change < 1.0 ) return true;
      do
      {
        types::t_unsigned i = rng.random(N-1);
        if ( xto_change > 1.0 and object.bitstring[i] < 0 )
          { object.bitstring[i] = 1; xto_change-=2; }
        else if ( xto_change < -1.0 and object.bitstring[i] > 0 )
          { object.bitstring[i] = -1; xto_change+=2; }

      } while ( xto_change < -1.0 or xto_change > 1.0 );
    }
    return true;
  }

}
