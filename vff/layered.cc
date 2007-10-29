//
//  Version: $Id$
//
#include <cstdlib>

#include "layered.h"
#include <opt/traits.h>

namespace Vff
{ 

  void Layered::create_template_strain()
  {
    // The first vector of the cell should indicate the direction of the
    // layering.
    u = structure.cell.get_column(0);
    types::t_real a = 1.0 / std::sqrt( atat::norm2(u) );
    u = a * u;
    template_strain.zero(); 

    // First, lets create an orthonormal vector to u
    atat::rVector3d a1 = Traits::Fuzzy<types::t_real>::equal( u(0), 0.0 ) ? 
                           ( Traits::Fuzzy<types::t_real>::equal( u(1), 0.0 ) ? 
                              atat::rVector3d(1, 0, 0):
                              atat::rVector3d(0, u(2), -u(3) )  ): 
                           atat::rVector3d( -u(2) -u(1), u(0), u(0) );
    a = ( 1.0 / std::sqrt( atat::norm2(a1) ) );
    a1 =  a * a1;

    // Then, lets create another... 
    atat::rVector3d a2;
    a2( 0 ) = u(1) * a1(2) - u(2) * a1(1);
    a2( 1 ) = u(2) * a1(0) - u(0) * a1(2);
    a2( 0 ) = u(0) * a1(1) - u(1) * a1(0);

    // Finally, we a transition matrix from a "u" basis to a cartesian basis.
    atat::rMatrix3d T; T.set_column(0, u); T.set_column(1, a1); T.set_column(2,a2);
    atat::rMatrix3d S; S.zero(); S(0,0) = 1;
    template_strain = T * S * (!T);
  }

  types::t_real Layered :: evaluate()
  {
    unpack_variables(strain);
    return t_Base::energy();
  }


  // initializes stuff before minimization
  bool Layered :: init()
  {
    create_template_strain();
    // sets up structure0, needed for fractional vs cartesian shit
    structure0 = structure;

    // Computes center of Mass
    // frozen (i.e. three components of the atomic positions )
    std::vector< Ising_CE::Atom > :: iterator i_atom =  structure0.atoms.begin();
    std::vector< Ising_CE::Atom > :: iterator i_atom_end =  structure0.atoms.end();
    center_of_mass = atat::rVector3d(0,0,0);
    for(; i_atom != i_atom_end; ++i_atom )
      center_of_mass += (!structure0.cell) * i_atom->pos;

    // Now counts the leftover degrees of freedom
    // There is already one assumed dof: relaxation of cell along the direction
    types::t_unsigned dof = 1;
    for( i_atom = structure0.atoms.begin(); 
         i_atom != i_atom_end; ++i_atom ) 
    {
      if ( not (i_atom->freeze & Ising_CE::Atom::FREEZE_X ) ) ++dof;
      if ( not (i_atom->freeze & Ising_CE::Atom::FREEZE_Y ) ) ++dof;
      if ( not (i_atom->freeze & Ising_CE::Atom::FREEZE_Z ) ) ++dof;
    }
    if ( not dof )
    {
      std::cerr << " Structure is frozen!! " << std::endl;
      std::cerr << " give me something to work with... " << std::endl;
      return false;
    }

    function::Base<> :: resize( dof );
    if ( not variables ) return false;

    strain.zero(); 
    strain(0,0) = 1.0;
    strain(1,1) = 1.0;
    strain(2,2) = 1.0;
    pack_variables(strain);
    
    return true;
  }

  // variables is expected to be of sufficient size!!
  void Layered :: pack_variables( const atat::rMatrix3d& _strain)
  {
    // finally, packs vff format into function::Base format
    iterator i_var = begin();
    *i_var = u * (_strain * u) - 1.0;
    ++i_var;

    std::vector< Ising_CE::Atom > :: const_iterator i_atom =  structure0.atoms.begin();
    std::vector< Ising_CE::Atom > :: const_iterator i_atom_end =  structure0.atoms.end();
    for(; i_atom != i_atom_end; ++i_atom )
    {
      if ( not (i_atom->freeze & Ising_CE::Atom::FREEZE_X ) )
        *i_var = i_atom->pos[0] * 0.5, ++i_var;
      if ( not (i_atom->freeze & Ising_CE::Atom::FREEZE_Y ) )
        *i_var = i_atom->pos[1] * 0.5, ++i_var;
      if ( not (i_atom->freeze & Ising_CE::Atom::FREEZE_Z ) )
        *i_var = i_atom->pos[2] * 0.5, ++i_var;
    }
  }


  // Unpacks opt::Function_Base::variables into Vff::Layered format
  void Layered :: unpack_variables(atat::rMatrix3d& strain)
  {
    std::vector<types::t_real> :: const_iterator i_x = variables->begin();
    std::vector<types::t_real> :: const_iterator i_x_end = variables->end();

    strain = (*i_x) * template_strain;
    strain(0,0) += 1.0;
    strain(1,1) += 1.0;
    strain(2,2) += 1.0;

    // compute resulting cell vectors
    structure.cell = strain * structure0.cell;
//   std::cout << " epsilon: " << *i_x << std::endl
//             << " template: " << std::endl << template_strain
//             << " strain: " << std::endl << strain
//             << " cell: " << std::endl << structure.cell << std::endl;

    // then computes positions
    std::vector<Ising_CE::Atom> :: const_iterator i_atom0 = structure0.atoms.begin();
    std::vector<Ising_CE::Atom> :: iterator i_atom = structure.atoms.begin();
    std::vector<Ising_CE::Atom> :: iterator i_atom_end = structure.atoms.end();
    atat::rVector3d com(0,0,0);
    atat::rMatrix3d cell_inv = !structure.cell;
    for(++i_x; i_atom != i_atom_end; ++i_atom, ++i_atom0 )
    {
      atat::rVector3d pos;
      pos[0] = ( i_atom0->freeze & Ising_CE::Atom::FREEZE_X ) ?
               i_atom0->pos[0] : 2.0 * (*i_x++);
      pos[1] = ( i_atom0->freeze & Ising_CE::Atom::FREEZE_Y ) ?
               i_atom0->pos[1] : 2.0 * (*i_x++);
      pos[2] = ( i_atom0->freeze & Ising_CE::Atom::FREEZE_Z ) ?
               i_atom0->pos[2] : 2.0 * (*i_x++);
      i_atom->pos = strain * pos;
      com -= cell_inv * i_atom->pos;
    }

    com += center_of_mass;
    com[0] /= (types::t_real) structure.atoms.size();
    com[1] /= (types::t_real) structure.atoms.size();
    com[2] /= (types::t_real) structure.atoms.size();
    if ( Traits::Fuzzy<types::t_real>::equal( atat::norm2( com ), 0 ) ) return;


    for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
      i_atom->pos += com; 
  }
} // namespace vff

