//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <algorithm>

#include <stdexcept>       // std::runtime_error

#include <physics/physics.h>
#include <opt/opt_frprmn.h>
#include <opt/debug.h>

#include "fortran.h"

//! Evaluation function for Minimizer::Frpr
extern "C" double fortranvff_for_frprmn( double* _x, double* _y)
    { return Minimizer::typical_frprfun<LaDa::Vff::Fortran>( _x, _y); }

namespace LaDa
{
  namespace Vff
  { 

    Fortran :: Fortran() : Functional(structure) 
    {
      Crystal::Structure::lattice = &lattice;

      lattice.cell.zero();
      lattice.cell(0,1) = 0.5;
      lattice.cell(0,2) = 0.5;
      lattice.cell(1,0) = 0.5;
      lattice.cell(1,2) = 0.5;
      lattice.cell(2,0) = 0.5;
      lattice.cell(2,1) = 0.5;
    
      Crystal::Lattice::t_Site site;
      site.site = 0; site.freeze = Crystal::Lattice::t_Site::FREEZE_NONE;
      site.pos = Eigen::Vector3d(0,0,0);
      lattice.sites.push_back(site);
      site.pos = Eigen::Vector3d(0.25,0.25,0.25);
      lattice.sites.push_back(site);
    }
    void Fortran::set_cell( types::t_real *_array )
    {
      // Watch it row vs column major. 
      types::t_real *i_var = _array;
      for(types::t_int i = 0; i < 3; ++i )
        for(types::t_int j = 0; j < 3; ++j, ++i_var )
          structure.cell( j, i ) = *i_var;
    }

    void Fortran :: set_atoms( types::t_int n, types::t_real *_position, char *_type )
    {
      std::vector< Crystal::StrAtom > atoms(n); 
      std::vector< Crystal::StrAtom > :: iterator i_atom = atoms.begin();
      std::vector< Crystal::StrAtom > :: iterator i_atom_end = atoms.end();

      types::t_real *i_pos = _position;
      char *i_char = _type;

      for(; i_atom != i_atom_end; ++i_atom, ++i_pos )
      {
        //! First determines position and lattice site
        i_atom->pos(0) = *i_pos;
        i_atom->pos(1) = *(++i_pos);
        i_atom->pos(2) = *(++i_pos);
        i_atom->site = lattice.get_atom_site_index(i_atom->pos);
        if( i_atom->site == -1 )
        {
          std::ostringstream sstr;
          sstr << __LINE__ << ", line: " << __LINE__ << "\n"
               << "Could not determine lattice-site of  " << i_atom->pos << std::endl;
        }
         

        //! Then extracts next symbol.
        std::string symbol;
        types::t_int n = Physics::Atomic::ExtractSymbol( i_char, symbol );
        __DOASSERT( n <= 0, "Error while extracting atomic symbol: " << symbol << "\n" )
        i_char += n;

        //! Then determines if it is a known type
        Crystal::Lattice::t_Site::t_Type &type = lattice.sites[ i_atom->site ].type;
        if( type.empty() )  type.push_back( symbol );
        else 
        {
          Crystal::Lattice::t_Site::t_Type :: const_iterator i_which;
          i_which = std::find( type.begin(), type.end(), symbol );
          if ( i_which == type.end() ) type.push_back( symbol );
        }

        // Finally sets the string atom
        i_atom->type = symbol;
        i_atom->freeze = Crystal::StrAtom::FREEZE_NONE;
      }

      // Finally, copies string atoms to structure
      for( i_atom = atoms.begin(); i_atom != i_atom_end; ++i_atom )
      {
        Crystal::Structure::t_Atom atom;
        lattice.convert_StrAtom_to_Atom( *i_atom, atom );
        structure.atoms.push_back( atom );
      }

      // creates an unitialized array of atomic functionals
      functionals.clear();
      functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(0,0),
                             structure, 0, 0) );
      if ( structure.lattice->get_nb_types(0) == 2 )
        functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(0,1),
                               structure, 0, 1) );
      functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(1,0),
                             structure, 1, 0) );
      if ( structure.lattice->get_nb_types(1) == 2 )
        functionals.push_back( Atomic_Functional( structure.lattice->get_atom_string(1,1),
                           structure, 1, 1) );

      //! Constructs the tree...
      initialize_centers();
    }

    void Fortran :: set_bond_parameters( char *_bond, const types::t_real _d0,
                                         const types::t_real *_alphas )
    {
      std::string A, B; 
      char *i_char = _bond;
      types::t_int n = Physics::Atomic::ExtractSymbol( i_char, A );
      __DOASSERT( n <= 0, "Error while extracting atomic symbol: " << A << "\n" )
      i_char += n;
      n = Physics::Atomic::ExtractSymbol( i_char, B );
      __DOASSERT( n <= 0, "Error while extracting atomic symbol: " << B << "\n" )
      i_char += n;

      
      types::t_int where[2];
      bond_indices( A, B, where );

      functionals[ where[0] ].add_bond( where[1], _d0, _alphas );
      functionals[ where[1]+structure.lattice->get_nb_types(0)].add_bond( where[0], 
                                                                          _d0, _alphas );
    }

    void Fortran :: set_angle_parameters( char *_angle, const types::t_real _gamma, 
                                          const types::t_real _sigma, const types::t_real *_betas )
    {
      std::string A, B, C; 
      char *i_char = _angle;
      types::t_int n = Physics::Atomic::ExtractSymbol( i_char, A );
      __DOASSERT( n <= 0, "Error while extracting atomic symbol: " << A << "\n" )
      i_char += n;
      n = Physics::Atomic::ExtractSymbol( i_char, B );
      __DOASSERT( n <= 0, "Error while extracting atomic symbol: " << B << "\n" )
      i_char += n;
      n = Physics::Atomic::ExtractSymbol( i_char, C );
      __DOASSERT( n <= 0, "Error while extracting atomic symbol: " << C << "\n" )
      i_char += n;

      types::t_int where[3];
      angle_indices( A, B, C, where );
      functionals[ where[1] ].add_angle( where[0], where[2], _gamma, _sigma, _betas );
    }
  }
} // namespace LaDa

extern "C" void FC_FUNC_(vff_create, VFF_CREATE)()
{
  if( LaDa::Vff :: fortran )
    __TRYDEBUGCODE( delete Vff :: fortran;,
                    "Could not delete Vff object.\n" )
  __TRYCODE( LaDa::Vff :: fortran = new LaDa::Vff::Fortran;,
             "Could not create Vff object.\n")
}
extern "C" void FC_FUNC_(vff_destroy, VFF_DESTROY)()
{
  if( LaDa::Vff :: fortran ) delete LaDa::Vff :: fortran;
  LaDa::Vff :: fortran = NULL;
}
extern "C" void FC_FUNC_(vff_print_structure, VFF_PRINT_STRUCTURE)()
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  LaDa::Vff :: fortran->print_structure();
}
extern "C" void FC_FUNC_(vff_print_lattice, VFF_PRINT_LATTICE)()
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  LaDa::Vff :: fortran->print_lattice();
}
extern "C" void FC_FUNC_(vff_scale, VFF_SCALE)( types::t_real* _scale )
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  LaDa::Vff :: fortran->set_scale( *_scale );
}
extern "C" void FC_FUNC_(vff_cell, VFF_CELL)( types::t_real* _mat)
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  LaDa::Vff :: fortran->set_cell( _mat );
}
extern "C" void FC_FUNC_(vff_atoms, VFF_ATOMS)( types::t_int *_n, 
                                                types::t_real* _pos, char *_type )
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  LaDa::Vff :: fortran->set_atoms( *_n, _pos, _type );
}
extern "C" void FC_FUNC_(vff_bond, VFF_BOND)( char *_bond, types::t_real *_d0,
                                              types::t_real *_alphas )
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  LaDa::Vff :: fortran->set_bond_parameters( _bond, *_d0, _alphas );
}
extern "C" void FC_FUNC_(vff_angle, VFF_ANGLE)( char *_angle, types::t_real *_gamma, 
                                                types::t_real *_sigma,
                                                types::t_real *_betas )
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  LaDa::Vff :: fortran->set_angle_parameters( _angle, *_gamma, *_sigma, _betas );
}
extern "C" void FC_FUNC_(vff_minimize, VFF_MINIMIZE)( types::t_real *_energy )
{
  __DOASSERT( not LaDa::Vff :: fortran, "Fortran object not initialized\n"; ) 
  typedef Minimizer::Frpr<LaDa::Vff::Fortran> t_Minimizer;
  t_Minimizer minimize( *LaDa::Vff :: fortran, fortranvff_for_frprmn );
  fortranvff_for_frprmn( (double*) LaDa::Vff::fortran, NULL );
  minimize.set_parameters( 4000, types::tolerance, 0.1, 0.01 );
  minimize();

  *_energy = LaDa::Vff::fortran->energy() / 16.0217733;
}
