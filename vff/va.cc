//
//  Version: $Id$
//
#include "va.h"

namespace Vff
{ 
  bool VirtualAtom :: Load( const TiXmlElement &_node )
  {
    // Load base
    if( not t_Base :: Load( _node ) ) return false;

    // Finds parent functional node
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
      if ( str.compare("vff" ) == 0 )
        break;
      parent = parent->NextSiblingElement("Functional");
    }
    if ( not parent )
    {
      std::cerr << "Could not find an <Functional type=\"vff\"> tag in input file" 
                << std::endl;
      return false;
    } 

    // Now checks for minimizer right in parent node
    if( minimizer.Load( *parent ) ) return true;

    // If it is not found within the functional, check for a minimizer node as
    // a sibling to parent;
    return minimizer.Load( *parent->Parent()->ToElement() );
  }

  bool VirtualAtom :: init()
  {
    va_vars.clear();
    va_vars.reserve( structure.atoms.size() );
    typedef Ising_CE :: Structure :: t_Atoms :: const_iterator t_ci;
    t_ci i_atom = structure.atoms.begin();
    t_ci i_atom_end = structure.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
      if( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) ) 
        va_vars.push_back( i_atom->type );

    return not va_vars.empty();
  }

  void VirtualAtom :: unpack_variables()
  {
    Ising_CE :: Structure :: t_Atoms :: iterator i_atom = structure.atoms.begin();
    Ising_CE :: Structure :: t_Atoms :: iterator i_atom_end = structure.atoms.end();
    t_Container :: const_iterator i_var = va_vars.begin();
    for(; i_atom != i_atom_end; ++i_atom )
    {
      if( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) continue;

      i_atom->type = *i_var > t_Type(0) ? 1.0: -1.0;
      ++i_var;
    }
  }

  VirtualAtom :: t_Type VirtualAtom :: evaluate()
  {
    unpack_variables();

    minimizer();
 
    structure.energy = energy();

    return structure.energy;
  }

  VirtualAtom :: t_Type VirtualAtom :: evaluate_one_gradient( types::t_unsigned _pos )
  {
    if( _pos > va_vars.size() )
      throw std::runtime_error( "Requesting out-of-range gradient.\n");

    std::vector< Atomic_Center > :: iterator i_center = centers.begin();
    std::vector< Atomic_Center > :: iterator i_center_end = centers.end();
    for(; _pos and i_center != i_center_end; ++i_center )
      if( not (i_center->Origin().freeze & Ising_CE::Structure::t_Atom::FREEZE_NONE) )
        --_pos;

    t_Type result = functionals[i_center->kind()].evaluate( *i_center );
    i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
    result -= functionals[i_center->kind()].evaluate( *i_center );
    result /= t_Type(2);
    i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);

    return result;
  }

  void VirtualAtom :: evaluate_gradient( t_Type * _grad )
  {
    t_Type* i_grad = _grad;
    std::vector< Atomic_Center > :: iterator i_center = centers.begin();
    std::vector< Atomic_Center > :: iterator i_center_end = centers.end();
    for(; i_center != i_center_end; ++i_center )
    {
      if( i_center->Origin().freeze & Ising_CE::Structure::t_Atom::FREEZE_NONE )
        continue;

      *i_grad = functionals[i_center->kind()].evaluate( *i_center );
      i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
      *i_grad -= functionals[i_center->kind()].evaluate( *i_center );
      *i_grad /= t_Type(2);
      i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);

      ++i_grad;
    } 

  }

  VirtualAtom :: t_Type VirtualAtom :: evaluate_with_gradient( t_Type * _grad )
  {
    t_Type result(0);
    t_Type* i_grad = _grad;
    std::vector< Atomic_Center > :: iterator i_center = centers.begin();
    std::vector< Atomic_Center > :: iterator i_center_end = centers.end();
    for(; i_center != i_center_end; ++i_center )
    {
      if( i_center->Origin().freeze & Ising_CE::Structure::t_Atom::FREEZE_NONE ) 
      {
        result += functionals[i_center->kind()].evaluate( *i_center );
        continue;
      }

      *i_grad = functionals[i_center->kind()].evaluate( *i_center );
      result += *i_grad;
      i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);
      *i_grad -= functionals[i_center->kind()].evaluate( *i_center );
      *i_grad /= t_Type(2);
      i_center->Origin().type = i_center->Origin().type > 0 ? t_Type(-1): t_Type(1);

      ++i_grad;
    } 

    return result;
  }

}
