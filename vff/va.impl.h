//
//  Version: $Id$
//
#ifndef _VFF_VA_IMPL_H_
#define _VFF_VA_IMPL_H_

namespace Vff
{

  template< class T_BASE >
  types::t_real PescanPosGrad :: position_gradient( types::t_unsigned _pos )
  {
    // Keep track of changes
    Ising_CE::Structure::t_Atoms copy_atoms = structure.atoms;
 
    // Flip atom and minimize strain
    Ising_CE::Structure::t_Atom& atom = structure.atoms[_pos];
    
    types::t_real result = vff.evaluate();

    // For periodicity reasons, we make the hypothesis that the cell-vectors
    // will not change with the atomic swap.
    atom.type = (atom.type > 0) ? -1.0: 1.0;
    types::t_unsigned cell_freeze = structure.freeze;
    structure.freeze = Ising_CE::Structure::FREEZE_ALL;
    vff.init();
    result -= vff.evaluate();
    structure.freeze = cell_freeze;
    atom.type = (atom.type > 0) ? -1.0: 1.0;

    // Finds out which atoms have changed sufficiently to warrant recomputation
    // And write em to file
    std::ostringstream stream;
    types::t_unsigned nb_pseudos=0; // nb of pseudos to put in perturbed potential

    Ising_CE::Structure :: const_iterator i_atom = atoms.begin();
    Ising_CE::Structure :: const_iterator i_atom_end = atoms.end();
    t_Centers :: const_iterator i_center = centers.begin();
    t_pseudos pseudos; pseudos.reserve(4);
    for(; i_atom != i_atom_end; ++i_atom, ++i_center )
    {
      //! Checks if this atom/center has moved appreciably
      typedef Ising_CE::Structure::t_Atoms::t_Atom t_atom;
      atat::rVector3d vec = ( i_out->pos - i_atom->pos ) * deriv_amplitude;
      if ( i_out->freeze & Ising_CE::Structure::t_Atom::FREEZE_X ) vec[0] = 0.0;
      if ( i_out->freeze & Ising_CE::Structure::t_Atom::FREEZE_Y ) vec[1] = 0.0;
      if ( i_out->freeze & Ising_CE::Structure::t_Atom::FREEZE_Z ) vec[2] = 0.0;
      if(     ( not Fuzzy<types::t_real> :: equal( vec[0], 0 ) )  
          and ( not Fuzzy<types::t_real> :: equal( vec[1], 0 ) )  
          and ( not Fuzzy<types::t_real> :: equal( vec[2], 0 ) )  ) continue;

      // then computes microscopic strain
      types::t_real msstrain = functionals[i_center->kind()].MicroStrain( *i_center,
                                                                          structure0 );
      
      // It has, hence it is printed out to the pseudo list.
      // First, we must find the type of empirical pseudo 
      find_escan_pseudo( i_center, pseudos );

      // now goes over found pseudos and creates output for both new position
      // (+derivative_amplitude) and old position (1-derivative_amplitude)
      t_pseudos::const_iterator i_pseudo = pseudos.begin();
      t_pseudos::const_iterator i_pseudo_end = pseudos.end();
      for( ; i_pseudo != i_pseudo_end; ++i_pseudo )
      {
        atat::rVector3d newpos = (!structure.cell) * i_center->Origin().pos;
        atat::rVector3d oldpos = (!structure.cell) * i_atom->pos;
        ++nb_pseudos;
               // old position
        stream << std::fixed    << std::setprecision(7)
               << std::setw(6)  << index << '0' << i_pseudo->first  // pseudo index
               << std::setw(12) << oldpos[0] // pseudo position
               << std::setw(12) << oldpos[1] 
               << std::setw(12) << oldpos[2] 
               << std::setw(18) << msstrain << " " // microscopic strain
               << std::setw(6) << std::setprecision(2)
                               <<   types::t_real( i_pseudo->second ) * 0.25  
                                  * (1.0 - deriv_amplitude)  // weight
               << std::setw(18) << std::setprecision(7) << pos[0] // pseudo position
               << std::setw(12) << pos[1] 
               << std::setw(12) << pos[2] << "\n"
               // new position
               << std::fixed    << std::setprecision(7)
               << std::setw(6)  << index << '0' << i_pseudo->first  // pseudo index
               << std::setw(12) << newpos[0] // pseudo position
               << std::setw(12) << newpos[1] 
               << std::setw(12) << newpos[2] 
               << std::setw(18) << msstrain << " " // microscopic strain
               << std::setw(6) << std::setprecision(2)
                               <<   types::t_real( i_pseudo->second ) * 0.25  
                                  * deriv_amplitude // weight
               << std::setw(18) << std::setprecision(7) << pos[0] // pseudo position
               << std::setw(12) << pos[1] 
               << std::setw(12) << pos[2] << "\n";
      }
    }

    // Opens file
    std::ofstream file( filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    // prints number of atoms
    file << nb_pseudos << "\n";
    // prints cell vectors in units of a0 and other
    // whatever nanopes other may be
    for( types::t_unsigned i = 0; i < 3; ++i )
      file << std::fixed << std::setprecision(7) 
           << std::setw(12) << structure.cell(0,i) * structure0.scale / Physics::a0("A")
           << std::setw(12) << structure.cell(1,i) * structure0.scale / Physics::a0("A")
           << std::setw(12) << structure.cell(2,i) * structure0.scale / Physics::a0("A")
           << std::setw(18) << structure.cell(0,i) 
           << std::setw(12) << structure.cell(1,i) 
           << std::setw(12) << structure.cell(2,i) << "\n";
    
    // print rest of file
    file << stream.str();
    file.flush();
    file.close();


    // Finally, copies structure back
    std::copy( copy_atoms.begin(), copy_atoms.end(),
               structure.atoms.begin() );

    return result/deriv_amplitude;
  }

  template< class T_BASE >
  types::t_real PescanPosGrad :: potential_gradient( types::t_unsigned _pos )
  {
    // Finds out which atoms have changed sufficiently to warrant recomputation
    // And write em to file
    std::ostringstream stream;
    types::t_unsigned nb_pseudos=0; // nb of pseudos to put in perturbed potential

    t_Centers :: const_iterator i_center = centers.begin() + _pos;
    types::t_real msstrain = functionals[i_center->kind()].MicroStrain( *i_center,
                                                                        structure0 );
    t_pseudos pseudos; pseudos.reserve(4);
      
    // first gets pseudo index
    Ising_CE::StrAtom stratom;
    structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[i_center->get_index()],
                                                stratom );
    types::t_unsigned index1 = Physics::Atomic::Z( stratom.type );
    // flips and gets other index


    atat::rVector3d centralpos = (!structure.cell) * i_center->Origin().pos;

    // then goes over bonds and prints pseudos
    Atomic_Center :: const_iterator i_bond = i_center->begin();
    Atomic_Center :: const_iterator i_bond_end = i_center->end();
    for(; i_bond != i_bond_end; ++i_bond )
    { 
      structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[i_bond->get_index()],
                                                  stratom );
      types::t_unsigned Z = Physics::Atomic::Z( stratom.type );

      ++nb_pseudos;
      
      // central "flipped" atom 
      stream << std::fixed    << std::setprecision(7)
             << std::setw(6)  << index << '0' << i_pseudo->first  // pseudo index
             << std::setw(12) << centralpos[0] // pseudo position
             << std::setw(12) << centralpos[1] 
             << std::setw(12) << centralpos[2] 
             << std::setw(18) << msstrain << " " // microscopic strain
             << std::setw(6) << std::setprecision(2)
                             <<   0.25 * (1.0 - deriv_amplitude)  // weight
             << std::setw(18) << std::setprecision(7) << pos[0] // pseudo position
             << std::setw(12) << pos[1] 
             << std::setw(12) << pos[2] << "\n";
    }
  }

  template< class T_BASE >
  void PescanPosGrad :: find_escan_pseudos( typename t_Centers::const_iterator &_i_center,
                                            t_pseudos _pseudos )
  {
      // first gets pseudo index
      Ising_CE::StrAtom stratom;
      structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[_i_center->get_index()],
                                                  stratom );
      types::t_unsigned index = Physics::Atomic::Z( stratom.type );

      
      // then goes over bonds and finds number of pseudos and their weights
      pseudos.clear();
      Atomic_Center :: const_iterator i_bond = _i_center->begin();
      Atomic_Center :: const_iterator i_bond_end = _i_center->end();
      for(; i_bond != i_bond_end; ++i_bond )
      { 
        structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[i_bond->get_index()],
                                                    stratom );
        types::t_unsigned Z = Physics::Atomic::Z( stratom.type );
        t_pseudos::iterator i_pseudo = _pseudos.begin();
        t_pseudos::iterator i_pseudo_end = _pseudos.end();
        for(; i_pseudo != i_pseudo_end; ++i_pseudo )
          if ( i_pseudo->first == Z ) break;

        if ( i_pseudo == i_pseudo_end ) _pseudos.push_back( t_pseudo( Z, 1 ) ); 
        else  ++(i_pseudo->second); 
      }
  }

} // namespace VFF

#endif // _VFF_VA_IMPL_H_
