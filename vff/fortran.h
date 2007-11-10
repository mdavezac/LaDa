//
//  Version: $Id$
//
#ifndef _FORTRAN_H_
#define _FORTRAN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>
#ifdef _MPI 
  #include <mpi/mpi_object.h>
#endif

#include <lamarck/structure.h>
#include <lamarck/lattice.h>

#include "functional.h"

extern "C" void vff_cell( types::t_real* );
extern "C" void vff_atoms( types::t_int *,  types::t_real*, char *_type );
extern "C" void vff_bonds( char *_bond, types::t_real *_d0,
                           types::t_real *_alphas ); 
extern "C" void vff_angles( char *_angle, types::t_real *_gamma, 
                            types::t_real *_sigma, types::t_real *_betas );
extern "C" void vff_create();
extern "C" void vff_destroy();
extern "C" void vff_minize( types::t_real *_energy );

namespace Vff
{
  class Fortran : public Functional
  {
    protected:
      Ising_CE::Structure structure;
      Ising_CE::Lattice lattice;

    public:
      Fortran();

      void set_cell( types::t_real *_array );
      void set_atoms( types::t_int n, types::t_real *_positions, char *_type );

      void set_bond_parameters( char *_bond, const types::t_real _d0,
                                const types::t_real *_alphas );
      void set_angle_parameters( char *_angle, const types::t_real _gamma, 
                                 const types::t_real _sigma, const types::t_real *_betas );
  };
 
  Fortran *fortran = NULL;
}

#endif
