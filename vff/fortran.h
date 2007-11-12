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

//! \ingroup Fortran
//! \brief Prints out vff structure
extern "C" void vff_print_structure_();
//! \ingroup Fortran
//! \brief Prints out vff lattice
extern "C" void vff_print_lattice_();
//! \ingroup Fortran
//! \brief Sets unit-cell of the vff structure.
//! \details Since fortran is the point, the array is a column-major 3x3
//!          matrix.
extern "C" void vff_cell_( types::t_real* );
//! \ingroup Fortran
//! \brief Sets the atomic position and the type of the vff structure.
//! \details This %function guesses the possible occupation of the lattice
//!          sites and creates the interaction tree of the structure. As a
//!          result, this function must be called before setting the
//!          interaction parameters.
//! \param _n Number of atoms_
//! \param _pos Array of types::t_real containing the cartesian coordinates of
//!             the atoms. In fortran, \a _pos is something like the following,
//!        \code
//  real*8 :: _pos( 3, _n )
//!        \endcode
//!             Whereas in C, the indices are reversed
//!        \code
//  types::t_real _pos[_n][3];
//!        \endcode
//!  \param _type Contains the types of each atom. The order should be the same
//!               as for pos. The string is expected to contain atomic symbols,
//!               separatated by blanks or '-' characters. In fortran, one can
//!               use an array of array of two characters.
//!         \code
//            character (len=2) :: At(_n)
//!         \endcode
extern "C" void vff_atoms_( types::t_int *_n,  types::t_real* _pos, char *_type );
//! \ingroup Fortran
//! \brief Sets the two-body interaction parameters of the bond described by \a _bond.
//! \param _bond chould an array characters with two atomic symbols separated
//!              by blanks or '-' characters.
//! \param _d0 is a real value representing the equilibrium bond-length.  
//! \param _alphas is an array of 5 real values representing the
//!                bond-stretching parameters. The first real value of the
//!                array should contain the quadratic bond-stretching
//!                parameter, the second real-value of the array should contain
//!                the cubic bond-stretching parameter, and so on.
extern "C" void vff_bond_( char *_bond, types::t_real *_d0,
                           types::t_real *_alphas ); 
//! \ingroup Fortran
//! \brief Sets the three body interactions of the three atoms described by \a _angle
//! \param _angle chould an array characters with three atomic symbols separated
//!              by blanks or '-' characters. The second atomic symbol
//!              represents the central atom of the interaction.
//! \param _gamma is a real value representing the (co?)sine of the equilibrium
//!               angle of the tetrahedra.
//! \param _sigma can't remember, Bond-angle parameter possibly.
//! \param _betas is an array of 5 real values representing the
//!                bond-bending parameters. The first real value of the
//!                array should contain the quadratic bond-bending
//!                parameter, the second real-value of the array should contain
//!                the cubic bond-bending parameter, and so on.
extern "C" void vff_angle_( char *_angle, types::t_real *_gamma, 
                            types::t_real *_sigma, types::t_real *_betas );
//! \ingroup Fortran
//! \brief Creates the vff object.
//! \details It is imperative to call this %function before all others.
extern "C" void vff_create_();
//! \ingroup Fortran
//! \brief Destroys the vff object.
//! \details This must be the last %function called.
extern "C" void vff_destroy_();
//! \ingroup Fortran
//! \brief Minimizes the structure and stores the \a _energy [eV].
extern "C" void vff_minimize_( types::t_real *_energy );
//! \ingroup Fortran
//! \brief Sets the scale of the cartesian units used in the structure and
//!        lattice.
extern "C" void vff_scale_( types::t_real* _scale );

namespace Vff
{
  //! \ingroup Fortran 
  //! \brief Interface object for fortran.
  //! \details This object is specialized for fcc lattices with two
  //!          atomic-sites (see Fortran::Fortran() ). An instance of this
  //!          object, pointed to by Vff::fortran, is created by vff_create_().
  //!          It is upon this instance that all the "C" functions of the
  //!          fortran interface will act.
  class Fortran : public Functional
  {
    protected:
      //! Ising_CE::Structure instance to minimize
      Ising_CE::Structure structure;
      //! Ising_CE::Lattice instance on which Fortran::structure exists.
      Ising_CE::Lattice lattice;

    public:
      //! Constructor
      Fortran();

      //! Prints the lattice to the standard output
      void print_lattice() { std::cout << lattice << std::endl; }
      //! Prints the structure to the standard output
      void print_structure() { std::cout << structure << std::endl; }
      //! Sets the scale of the structure and the lattice
      void set_scale( types::t_real _scale )
        { structure.scale = lattice.scale = _scale;  }
      //! \brief Sets the cell of the structure.
      //! \details Since fortran is the point, the array is a column-major 3x3
      //!          matrix.
      void set_cell( types::t_real *_array );
      //! \brief Sets the atomic positions.
      //! \see vff_atom_()
      void set_atoms( types::t_int n, types::t_real *_positions, char *_type );

      //! \brief sets the parameters of bond interactions.
      //! \see vff_bond()
      void set_bond_parameters( char *_bond, const types::t_real _d0,
                                const types::t_real *_alphas );
      
      //! \brief sets the parameters of angle and bond-angle interactions.
      //! \see vff_angle()
      void set_angle_parameters( char *_angle, const types::t_real _gamma, 
                                 const types::t_real _sigma, const types::t_real *_betas );
  };
 
  //! Global instance on which the fortran interface will act.
  Fortran *fortran = NULL;
}

#endif
