//
//  Version: $Id$
//

#ifndef _LAYERED_IMPL_H_
#define _LAYERED_IMPL_H_


#include <print/stdout.h>
#include <print/manip.h>
#include <lamarck/structure.h>
#include <lamarck/atom.h>

namespace Layered
{

  template <types::t_unsigned _D>
  inline Concentration<_D> :: Concentration () : N(0), single_c(false) 
  {
    for(unigned i=0; i < _D; ++i )
     { x0[i] = 0; x[i] = 0; Nfreeze[i] = 0; }
  }
  template <types::t_unsigned _D>
  inline Concentration<_D> :: Concentration   ( const Concentration &_conc)
                                            : N(_conc.N), single_c(_conc.single_c) 
  {
    for(unigned i=0; i < _D; ++i )
      { x0[i] = _conc.x0[i]; x[i] = _conc.x[i]; Nfreeze[i] = _conc.Nfreeze[i]; }
  }


  template <types::t_unsigned _D>
  void Concentration<_D> :: operator()( Ising_CE::Structure &_str )
  {
    if ( _str.atoms.size() % _D != 0 )
      throw std::runtime_error( "Number of atoms is not a multiple of _D in
                                 Layered::Concentration\n");

    if ( not single_c ) return;

    types::t_complex  *hold = new types::t_complex[ N ];
    if ( not hold )
    {
      std::cerr << " Could not allocate memory in set_concentration!! " << std::endl; 
      exit(0);
    }

    // creates individual with unnormalized occupation
    types::t_complex  *i_hold = hold;
    Fourier<_D>( _str.atoms.begin(), _str.atoms.end(),
                 _str.k_vecs.begin(), _str.k_vecs.end(),
                 i_hold );

    Ising_CE::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Ising_CE::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int concx = 0;
    i_hold = hold;
    if( not single_c )
    {
      for (; i_atom != i_atom_end; i_atom += _D, ++i_hold)
        if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        {
          if ( std::abs( std::real(*i_hold) ) < types::tolerance )
                i_atom->type = rng.flip() ? 1.0: -1.0;
          else  i_atom->type = std::real(*i_hold) > 0.0 ? 1.0: -1.0;
          (i_atom+1)->type = i_atom->type > 0.0 ? 1.0: -1.0;
        }
      return;
    }

    for (; i_atom != i_atom_end; i_atom+=_D, ++i_hold)
    {
      if ( not ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
        i_atom->type = std::real(*i_hold);
      ( i_atom->type > 0 ) ? ++concx : --concx;
    }

    // then normalize it while setting correct concentrations
    normalize( _str, (types::t_real) concx - ( (types::t_real) N ) * x0 );
    delete[] hold;
    return;
  }



  // Takes an "unphysical" individual and set normalizes its sites _sites to +/-1,
  // after flipping the _tochange spins closest to zero.
  // ie sets the concentration
  template <types::t_unsigned _D>
  void Concentration<_D> :: normalize( Ising_CE::Structure &_str, types::t_real _tochange ) 
  {
    Ising_CE::Structure::t_Atoms::iterator i_end = _str.atoms.end();
    Ising_CE::Structure::t_Atoms::iterator i_which;
    Ising_CE::Structure::t_Atoms::iterator i_atom;
    while( _tochange < -1.0 or _tochange > 1.0 )
    {
      // first finds atom with type closest to zero from +1 or -1 side,
      // depending on _tochange
      i_atom = _str.atoms.begin();
      i_which = i_end;
      types::t_real minr = 0.0;
      for(; i_atom != i_end; i_atom += _D )
      {
        if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) continue;
        if ( _tochange > 0 )
        {
          if ( i_atom->type < 0 )                    continue;
          if ( minr != 0.0 and i_atom->type > minr ) continue;
        }
        else // ( _tochange < 0 )
        {
          if ( i_atom->type > 0 )                    continue;
          if ( minr != 0.0 and i_atom->type < minr ) continue;
        }

        i_which = i_atom;
        minr = i_atom->type;
      }
      if ( i_which == i_end )
        throw std::runtime_error( "Error while normalizing x constituents\n" );

      i_which->type = ( _tochange > 0 ) ? -1.0: 1.0;
      _tochange += ( _tochange > 0 ) ? -2: 2;
    }

    // finally, normalizes _str
    i_atom = _str.atoms.begin();
    for(unsigned i=0; i_atom != i_end; ++i_atom, ++i )
    {
      i_atom->type = ( i_atom->type > 0 ) ? 1.0: -1.0;
      i_which = i_atom + 1;
      for(types::t_unsigned i = 0; i < _D; ++i )
        i_which->type = i_atom->type;
    }
  }








} // namespace Layered


#endif
