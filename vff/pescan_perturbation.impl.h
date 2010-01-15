//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "pescan_perturbation.h"

#include <physics/physics.h>

namespace LaDa
{
  namespace Vff
  { 
    template< class T_VFFBASE >
    types::t_real PescanPerturbations<T_VFFBASE> :: stress( types::t_unsigned _pos )
    {
      // Keep track of changes
      t_Atoms copy_atoms = t_VA::structure.atoms;
      // Atom which will be changed
      t_Atom& atom = t_VA::structure.atoms[_pos];

      // Minimize strain and evaluate 
      types::t_real result = evaluate();

      // Original MicroStrain
      std::vector< types::t_real > microstrain;
      microstrain.reserve( t_VA::structure.atoms.size() );
      t_Atoms :: const_iterator i_atom = t_VA::structure.atoms.begin();
      t_Atoms :: const_iterator i_atom_end = t_VA::structure.atoms.end();
      t_Center& center = centers[_pos];
      for(; i_atom != i_atom_end; ++i_atom )
        microstrain.push_back( functionals[center.kind()]
                                          .MicroStrain( center, structure0 ) );

      // For periodicity reasons, we make the hypothesis that the cell-vectors
      // will not change with the atomic swap.
      types::t_unsigned cell_freeze = t_VA::structure.freeze;
      t_VA::structure.freeze = Crystal::Structure::FREEZE_ALL;

      // Flip atom and minimize strain
      atom.type = (atom.type > 0) ? -1.0: 1.0;
      t_VA::init(); 
      result -= evaluate();
      // undo changes to structure now
      t_VA::structure.freeze = cell_freeze;
      atom.type = (atom.type > 0) ? -1.0: 1.0;

      // Finds out which atoms have changed sufficiently to warrant recomputation
      // And write em to file
      std::ostringstream stream;
      types::t_unsigned nb_pseudos=0; // nb of pseudos to put in perturbed potential

      i_atom = t_VA::structure.atoms.begin();
      typename t_Centers :: const_iterator i_center = centers.begin();
      std::vector<types::t_real> :: const_iterator i_ms = microstrain.begin();
      t_pseudos pseudos; pseudos.reserve(4);
      for(; i_atom != i_atom_end; ++i_atom, ++i_center, ++i_ms )
      {
        //! Checks if this atom/center has moved appreciably
        typedef t_Atom t_atom;
        Eigen::Vector3d vec = ( i_center->Origin().pos - i_atom->pos ) * deriv_amplitude;
        // then computes microscopic strain
        types::t_real msstrain = functionals[i_center->kind()].MicroStrain( *i_center,
                                                                            structure0 );
        // then goes over bonds and finds number of pseudos and their weights
        if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_X ) vec[0] = 0.0;
        if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_Y ) vec[1] = 0.0;
        if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_Z ) vec[2] = 0.0;
        if(     ( not Fuzzy::eq( vec[0], 0e0 ) )  
            and ( not Fuzzy::eq( vec[1], 0e0 ) )  
            and ( not Fuzzy::eq( vec[2], 0e0 ) )  ) continue;

        // It has, hence it is printed out to the pseudo list.
        // First, we must find the type of empirical pseudo 
        find_escan_pseudos( *i_center, pseudos );

        // gets pseudo Z
        Crystal::StrAtom stratom;
        t_VA::structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[i_center->get_index()],
                                                          stratom );
        types::t_unsigned Z = Physics::Atomic::Z( stratom.type );

        
        // now goes over found pseudos and creates output for both new position
        // (+derivative_amplitude) and old position (1-derivative_amplitude)
        t_pseudos::const_iterator i_pseudo = pseudos.begin();
        t_pseudos::const_iterator i_pseudo_end = pseudos.end();
        for( ; i_pseudo != i_pseudo_end; ++i_pseudo )
        {
          Eigen::Vector3d newpos = (!t_VA::structure.cell) * i_center->Origin().pos;
          Eigen::Vector3d oldpos = (!t_VA::structure.cell) * i_atom->pos;
          nb_pseudos += 2;
                 // old position and microstrain
          stream << std::fixed    << std::setprecision(7)
                 << std::setw(6)  << Z << '0' << i_pseudo->first  // pseudo index
                 << std::setw(12) << oldpos[0] // pseudo position
                 << std::setw(12) << oldpos[1] 
                 << std::setw(12) << oldpos[2] 
                 << std::setw(18) << *i_ms << " " // microscopic strain
                 << std::setw(6) << std::setprecision(2)
                                 <<   types::t_real( i_pseudo->second ) * 0.25  
                                    * (1.0 - deriv_amplitude)  // weight
                 << std::setw(18) << std::setprecision(7) << oldpos[0] // pseudo position
                 << std::setw(12) << oldpos[1] 
                 << std::setw(12) << oldpos[2] << "\n"
                 // new position and microstrain
                 << std::fixed    << std::setprecision(7)
                 << std::setw(6)  << Z << '0' << i_pseudo->first  // pseudo index
                 << std::setw(12) << newpos[0] // pseudo position
                 << std::setw(12) << newpos[1] 
                 << std::setw(12) << newpos[2] 
                 << std::setw(18) << msstrain << " " // microscopic strain
                 << std::setw(6) << std::setprecision(2)
                                 <<   types::t_real( i_pseudo->second ) * 0.25  
                                    * deriv_amplitude // weight
                 << std::setw(18) << std::setprecision(7) << newpos[0] // pseudo position
                 << std::setw(12) << newpos[1] 
                 << std::setw(12) << newpos[2] << "\n";
        }
      }

      // Opens file
      std::ofstream file( filename.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 
      // prints number of atoms
      file << nb_pseudos << "\n";
      // prints cell vectors in units of a0 and other
      // whatever nanopes other may be
      for( types::t_unsigned i = 0; i < 3; ++i )
        file << std::fixed << std::setprecision(7) 
             << std::setw(12) <<   t_VA::structure.cell(0,i)
                                 * structure0.scale / Physics::a0("A")
             << std::setw(12) <<   t_VA::structure.cell(1,i)
                                 * structure0.scale / Physics::a0("A")
             << std::setw(12) <<   t_VA::structure.cell(2,i)
                                 * structure0.scale / Physics::a0("A")
             << std::setw(18) << t_VA::structure.cell(0,i) 
             << std::setw(12) << t_VA::structure.cell(1,i) 
             << std::setw(12) << t_VA::structure.cell(2,i) << "\n";
      
      // print rest of file
      file << stream.str();
      file.flush();
      file.close();


      // Finally, copies structure back
      std::copy( copy_atoms.begin(), copy_atoms.end(),
                 t_VA::structure.atoms.begin() );

      return result/deriv_amplitude;
    }

    template< class T_VFFBASE >
    void PescanPerturbations<T_VFFBASE> :: chemical( types::t_unsigned _pos )
    {
      // Only the tetrahedron around pos changes
      std::ostringstream stream;
      types::t_unsigned nb_pseudos=0; // nb of pseudos to put in perturbed potential

      t_Center &center = centers[_pos];
      types::t_real msstrain = functionals[center.kind()].MicroStrain( center,
                                                                       structure0 );
      t_pseudos pseudos; pseudos.reserve(4);
        
      // first gets pseudo index
      Crystal::StrAtom stratom;
      t_VA::structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[center.get_index()],
                                                        stratom );
      types::t_unsigned index1 = Physics::Atomic::Z( stratom.type );
      // flips and gets other index
      center.Origin().type = ( center.Origin().type > 0 ? -1.0: 1.0 );
      t_VA::structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[center.get_index()],
                                                        stratom );
      types::t_unsigned index2 = Physics::Atomic::Z( stratom.type );
      // flips back before all hell breaks loose
      center.Origin().type = ( center.Origin().type > 0 ? -1.0: 1.0 );



      //! central position
      Eigen::Vector3d centralpos = (!t_VA::structure.cell) * center.Origin().pos;

      // then goes over bonds and prints pseudos
      AtomicCenter :: const_iterator i_bond = center.begin();
      AtomicCenter :: const_iterator i_bond_end = center.end();
      for(; i_bond != i_bond_end; ++i_bond )
      { 
        t_VA::structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[i_bond->get_index()],
                                                          stratom );
        types::t_unsigned Z = Physics::Atomic::Z( stratom.type );

        Eigen::Vector3d bondpos = (!t_VA::structure.cell) * i_bond->Origin().pos;

        nb_pseudos += 4;
        // central "original" atom 
        stream << std::fixed    << std::setprecision(7)
               << std::setw(6)  << index1 << '0' << Z  // pseudo index
               << std::setw(12) << centralpos[0] // pseudo position
               << std::setw(12) << centralpos[1] 
               << std::setw(12) << centralpos[2] 
               << std::setw(18) << msstrain << " " // microscopic strain
               << std::setw(6) << std::setprecision(2)
                               <<   0.25 * (1.0 - deriv_amplitude)  // weight
               << std::setw(18) << std::setprecision(7)  // pseudo position
                                << centralpos[0]
               << std::setw(12) << centralpos[1] 
               << std::setw(12) << centralpos[2] << "\n"
        //  central "flipped" atom
               << std::fixed    << std::setprecision(7)
               << std::setw(6)  << index2 << '0' << Z  // pseudo index
               << std::setw(12) << centralpos[0] // pseudo position
               << std::setw(12) << centralpos[1] 
               << std::setw(12) << centralpos[2] 
               << std::setw(18) << msstrain << " " // microscopic strain
               << std::setw(6) << std::setprecision(2)
                               <<   0.25 * deriv_amplitude  // weight
               << std::setw(18) << std::setprecision(7)  // pseudo position
                                << centralpos[0]
               << std::setw(12) << centralpos[1] 
               << std::setw(12) << centralpos[2] << "\n"
        // bond "original" atom 
               << std::fixed    << std::setprecision(7)
               << std::setw(6)  << Z << '0' << index1  // pseudo index
               << std::setw(12) << bondpos[0] // pseudo position
               << std::setw(12) << bondpos[1] 
               << std::setw(12) << bondpos[2] 
               << std::setw(18) << msstrain << " " // microscopic strain
               << std::setw(6) << std::setprecision(2)
                               <<   0.25 * (1.0 - deriv_amplitude)  // weight
               << std::setw(18) << std::setprecision(7)  // pseudo position
                                << bondpos[0]
               << std::setw(12) << bondpos[1] 
               << std::setw(12) << bondpos[2] << "\n"
        //  bond "flipped" atom
               << std::fixed    << std::setprecision(7)
               << std::setw(6)  << Z << '0' << index2  // pseudo index
               << std::setw(12) << bondpos[0] // pseudo position
               << std::setw(12) << bondpos[1] 
               << std::setw(12) << bondpos[2] 
               << std::setw(18) << msstrain << " " // microscopic strain
               << std::setw(6) << std::setprecision(2)
                               <<   0.25 * deriv_amplitude  // weight
               << std::setw(18) << std::setprecision(7)  // pseudo position
                                << bondpos[0]
               << std::setw(12) << bondpos[1] 
               << std::setw(12) << bondpos[2] << "\n";
      }

      // Opens file
      std::ofstream file( filename.string().c_str(), std::ios_base::out|std::ios_base::trunc ); 
      // prints number of atoms
      file << nb_pseudos << "\n";
      // prints cell vectors in units of a0 and other
      // whatever nanopes other may be
      for( types::t_unsigned i = 0; i < 3; ++i )
        file << std::fixed << std::setprecision(7) 
             << std::setw(12) <<   t_VA::structure.cell(0,i)
                                 * structure0.scale / Physics::a0("A")
             << std::setw(12) <<   t_VA::structure.cell(1,i)
                                 * structure0.scale / Physics::a0("A")
             << std::setw(12) <<   t_VA::structure.cell(2,i)
                                 * structure0.scale / Physics::a0("A")
             << std::setw(18) << t_VA::structure.cell(0,i) 
             << std::setw(12) << t_VA::structure.cell(1,i) 
             << std::setw(12) << t_VA::structure.cell(2,i) << "\n";
      
      // print rest of file
      file << stream.str();
      file.flush();
      file.close();
    }

    template<class T_VFFBASE>
    void PescanPerturbations<T_VFFBASE> :: find_escan_pseudos( const t_Center &_center,
                                                               t_pseudos _pseudos )
    {
      _pseudos.clear();
      typename t_Center :: const_iterator i_bond = _center.begin();
      typename t_Center :: const_iterator i_bond_end = _center.end();
      for(; i_bond != i_bond_end; ++i_bond )
      { 
        Crystal::StrAtom stratom;
        t_VA::structure.lattice->convert_Atom_to_StrAtom( structure0.atoms[i_bond->get_index()],
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

  }
} // namespace LaDa
