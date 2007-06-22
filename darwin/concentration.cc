#include "concentration.h"

types::t_real set_concentration( Ising_CE::Structure &_str,
                                 types::t_real _target )
{
  types::t_unsigned N = (types::t_int) _str.atoms.size();
  types::t_complex  *hold = new types::t_complex[ N ];
  if ( not hold )
  {
    std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
    exit(0);
  }
  types::t_complex  *i_hold = hold;
  Ising_CE::fourrier_to_rspace( _str.atoms.begin(), _str.atoms.end(),
                                _str.k_vecs.begin(), _str.k_vecs.end(),
                                i_hold );

  Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
  Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
  types :: t_int result = 0; 
  i_hold = hold;
  for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
  {
    if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
      ( i_atom->type > 0 ) ? ++result : --result;
    else 
    {
      ( std::real(*i_hold) > 0 ) ? ++result : --result;
      i_atom->type = ( std::real(*i_hold) > 0 ) ? 1.0: -1.0;
    }
  }

  if ( _target == -2.0 )
  {
    delete[] hold;
    return (double) result  / (double) N;
  }

  types::t_int to_change = static_cast<types::t_int>( (double) N * _target ) - result;
  if (  not to_change ) 
  {
    delete[] hold;
    return _target;
  }

  i_atom = _str.atoms.begin();
  i_hold = hold;
  for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
      *i_hold = to_change;
    else if ( to_change > 0 and std::real( *i_hold ) > 0 )
      *i_hold = to_change;
    else if ( to_change < 0 and std::real( *i_hold ) < 0 )
      *i_hold = to_change;

  i_atom = _str.atoms.begin();
  types::t_int which;
  do
  {
    i_hold = hold;
    which = 0;
    types::t_real maxr = 1.0;
    if ( to_change > 0 )
    {
      for( types::t_unsigned i = 0; i < N; ++i, ++i_hold ) 
        if (     ( maxr == 1.0 and std::real( *i_hold ) < 0 )
             or  (  std::real( *i_hold ) < 0 and maxr <  std::real( *i_hold ) ) )
        {
          maxr = std::real( *i_hold );
          which = i;
        }
      ( i_atom + which )->type = 1.0;
      *( hold + which ) = to_change;
      result+=2; to_change-=2;
      continue;
    }
    maxr = -1.0;
    for( types::t_unsigned i = 0; i < N; ++i, ++i_hold ) 
      if (     ( maxr == -1.0 and std::real( *i_hold ) > 0 )
           or  (  std::real( *i_hold ) > 0 and maxr >  std::real( *i_hold ) ) )
      {
        maxr = std::real( *i_hold );
        which = i;
      }
    ( i_atom + which )->type = -1.0;
    *( hold + which ) = to_change;
    result-=2; to_change+=2;

  } while (to_change); 
  
  delete[] hold;
  return (double) result  / (double) N;
}


  bool X_vs_y :: Load ( const TiXmlElement &_node )
  {
    const TiXmlElement *parent;
    std::string str;

    // This whole section tries to find a <Functional type="vff"> tag
    // in _node or its child
    str = _node.Value();
    if ( str.compare("Functional" ) != 0 )
      parent = _node.FirstChildElement("Functional");
    else
      parent = &_node;
    while (parent)
    {
      str = "";
      if ( parent->Attribute( "type" )  )
        str = parent->Attribute("type");
      if ( str.compare("Concentration" ) == 0 )
        break;
      parent = parent->NextSiblingElement("Functional");
    }
    if ( not parent )
    {
      std::cerr << "Could not find an <Functional type=\"Concentration\"> tag in input file" 
                << std::endl;
      return false;
    } 

    if (    ( not parent->Attribute("a") )
         or ( not parent->Attribute("b") )
         or ( not parent->Attribute("c") ) )
    {
      std::cerr << "Concentration functional is missing some input "<<std::cerr;
      return false;
    }
    parent->Attribute("a", &a);
    parent->Attribute("b", &b);
    parent->Attribute("c", &c);

    if (    ( b*b - 4*a*c < 0.0 )
         or ( b*b - 4*(a-1.0)*c < 0.0 ) )
    {
      std::cerr << "Equation incorrect on input... should be y = a +b*x + c*x*x "<<std::cerr;
      std::cerr << " with x,y in [0,1] "<<std::cerr;
    }

    // changing to x,y in [-1,1]
    a += b *0.5 + c*0.25 -0.5;
    b = b *0.5 + c*0.5;
    c *= 0.25;
    return true;
  }

#ifdef _MPI
namespace mpi
{
  template<>
  bool mpi::BroadCast::serialize<X_vs_y>( X_vs_y & _xy )
  {
    if( not serialize( _xy.a ) ) return false;
    if( not serialize( _xy.b ) ) return false;
    if( not serialize( _xy.c ) ) return false;

    return true;
  }
}
#endif

