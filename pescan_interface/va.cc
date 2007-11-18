//
//  Version: $Id$
//
#include "va.h"

namespace Pescan
{ 
  bool VirtualAtom :: Load( const TiXmlElement &_node )
  {
    // Load base
    return     pescan.Load( _node ) 
           and vff.Load( _node );
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

  void VirtualAtom :: position_grad( types::t_unsigned _pos )
  {
    // Keep track of changes
    Ising_CE::Structure copy_structure = structure;
    Bands copy_bands = bands;
 
    // Flip atom and minimize strain
    Ising_CE::Structure::t_Atom& atom = structure.atoms[_pos];
    
    atom.type = (atom.type > 0) ? -1.0: 1.0;
    vff.init();
    vff.evaluate();
    atom.type = (atom.type > 0) ? -1.0: 1.0;

    // Rewrites structure with changed position
    Ising_CE::Structure::t_Atoms::const_iterator i_atom = copy_structure.begin();
    Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = copy_structure.end();
    Ising_CE::Structure::t_Atoms::iterator i_out = copy_structure.begin();
    for(; i_atom != i_atom_end; ++i_atom, ++i_out )
    {
      typedef Ising_CE::Structure::t_Atoms::t_Atom t_atom;
      atat::rVector3d vec = ( i_out->pos - i_atom->pos ) * deriv_amplitude;
      if(     ( not ( i_atom->freeze & t_Atom :: FREEZE_X ) )
          and ( not Fuzzy<types::t_real> :: equal( vec[0], 0 ) ) )
        i_out->pos[0] =  vec[0] + i_atom->out[0];
      if(     ( not ( i_atom->freeze & t_Atom :: FREEZE_Y ) )
          and ( not Fuzzy<types::t_real> :: equal( vec[1], 0 ) ) )
        i_out->pos[1] =  vec[1] + i_atom->out[1];
      if(     ( not ( i_atom->freeze & t_Atom :: FREEZE_Z ) )
          and ( not Fuzzy<types::t_real> :: equal( vec[2], 0 ) ) )
        i_out->pos[2] =  vec[2] + i_atom->out[2];
    }

    // Now computes bandgap

  }

}
