//
//  Version: $Id$
//
#include "va.h"

namespace Vff
{ 
  bool VirtualAtom :: init()
  {
    va_vars.clear();
    va_vars.reserve( structure.atoms.size() );

    Ising_CE :: Structure :: t_Atoms :: const_iterator i_atom = structure.begin();
    Ising_CE :: Structure :: t_Atoms :: const_iterator i_atom_end = structure.end();
    for(; i_atom != i_atom_end; ++i_atom )
      if( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) ) 
        va_vars.push_back( i_atom->type );

    return not va_vars.empty();
  }

  void VirtualAtom :: unpack_variables()
  {
    Ising_CE :: Structure :: t_Atoms :: iterator i_atom = structure.begin();
    Ising_CE :: Structure :: t_Atoms :: iterator i_atom_end = structure.end();
    t_Container :: const_iterator i_var = va_vars.begin();
    for(; i_atom != i_atom_end; ++i_atom )
    {
      if( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) continue;

      i_atom->type = *i_var > t_Type(0) ? 1.0: -1.0;
      ++i_var
    }
  }

  VirtualAtom :: t_Type VirtualAtom :: evaluate()
  {
    unpack_variables();

    minimize();
 
    structure.energy = energy();

    return structure.energy;
  }

  VirtualAtom :: t_Type VirtualAtom :: evaluate_one_gradient( types::t_unsigned _pos )
  {
    if( _pos > va_vars.size() ) throw std::runtime_error( "Requesting out-of-range gradient.\n");

    std::vector< Atomic_Center > :: iterator i_center = centers.begin();
    std::vector< Atomic_Center > :: iterator i_center_end = centers.end();
    for(; _pos and i_center != i_center_end; ++i_center )
      if( not ( i_center->origin->freeze & Ising_CE::Structure::t_Atom::FREEZE_NONE ) ) --_pos;

    t_Type result = i_center->evaluate();
    i_center->origin->type = i_center->origin->type > 0 ? t_Type(-1): t_Type(1);
    result -= i_center->evaluate();
    result /= t_Type(2);
    i_center->origin->type = i_center->origin->type > 0 ? t_Type(-1): t_Type(1);

    return result;
  }

  void VirtualAtom :: evaluate_gradient( t_Type * _grad )
  {
    t_Type* i_grad = _grad;
    std::vector< Atomic_Center > :: iterator i_center = centers.begin();
    std::vector< Atomic_Center > :: iterator i_center_end = centers.end();
    for(; _pos and i_center != i_center_end; ++i_center )
    {
      if( i_center->origin->freeze & Ising_CE::Structure::t_Atom::FREEZE_NONE ) continue;

      *i_grad = functionals[i_center->kind()].evaluate( *i_center );
      i_center->origin->type = i_center->origin->type > 0 ? t_Type(-1): t_Type(1);
      *i_grad -= functionals[i_center->kind()].evaluate( *i_center );
      *i_grad /= t_Type(2);
      i_center->origin->type = i_center->origin->type > 0 ? t_Type(-1): t_Type(1);

      ++i_grad;
    } 

  }

  VirtualAtom :: t_Type VirtualAtom :: evaluate_with_gradient( t_Type * _grad )
  {
    t_Type result(0);
    t_Type* i_grad = _grad;
    std::vector< Atomic_Center > :: iterator i_center = centers.begin();
    std::vector< Atomic_Center > :: iterator i_center_end = centers.end();
    for(; _pos and i_center != i_center_end; ++i_center )
    {
      if( i_center->origin->freeze & Ising_CE::Structure::t_Atom::FREEZE_NONE ) 
      {
        result += functionals[i_center->kind()].evaluate( *i_center );
        continue;
      }

      *i_grad = functionals[i_center->kind()].evaluate( *i_center );
      result += *i_grad;
      i_center->origin->type = i_center->origin->type > 0 ? t_Type(-1): t_Type(1);
      *i_grad -= functionals[i_center->kind()].evaluate( *i_center );
      *i_grad /= t_Type(2);
      i_center->origin->type = i_center->origin->type > 0 ? t_Type(-1): t_Type(1);

      ++i_grad;
    } 

    return result;
  }

}

#endif
