//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdlib>

#include <algorithm>
#include <functional>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/filesystem/operations.hpp>
#ifdef _MPI
# include <boost/mpi/collectives/all_reduce.hpp>
#endif

#include <physics/physics.h>
#include <opt/atat.h>
#include <opt/debug.h>
#include <opt/atat.h>
#include <opt/tinyxml.h>

#include "functional.h"
  
namespace LaDa
{
  namespace Vff
  { 
    Functional :: t_Return Functional :: operator()( const t_Arg& _arg ) const
    {
      atat::rMatrix3d strain;
      unpack_variables( _arg, strain );
      return Vff::energy();
    }

    // Unpacks opt::Function_Base::variables into Vff::Functional format
    void Functional :: unpack_variables(const t_Arg& _arg, atat::rMatrix3d& _strain) const
    {
      t_Arg :: const_iterator i_x = _arg.begin();

      _strain(0,0) = ( structure.freeze & Crystal::Structure::FREEZE_XX ) ?
                       types::t_real(1.0) : (*i_x++);
      _strain(1,1) = ( structure.freeze & Crystal::Structure::FREEZE_YY ) ?
                       types::t_real(1.0) : (*i_x++);
      _strain(2,2) = ( structure.freeze & Crystal::Structure::FREEZE_ZZ ) ?
                       types::t_real(1.0) : (*i_x++);
      _strain(0,1) = _strain (1,0) = (structure.freeze & Crystal::Structure::FREEZE_XY) ?
                                       types::t_real(0.0) : (*i_x++);
      _strain(0,2) = _strain (2,0) = (structure.freeze & Crystal::Structure::FREEZE_XZ) ?
                                       types::t_real(0.0) : (*i_x++);
      _strain(2,1) = _strain (1,2) = (structure.freeze & Crystal::Structure::FREEZE_YZ) ?
                                       types::t_real(0.0) : (*i_x++);

      // compute resulting cell vectors
      structure.cell = _strain * structure0.cell;
      unpack_positions( i_x, _strain );
    }

    void Functional :: unpack_positions( t_Arg :: const_iterator &_i_x, 
                                         atat::rMatrix3d& _strain ) const
    {
      // then computes positions
      t_Atoms :: iterator i_atom = structure.atoms.begin();
      t_Atoms :: iterator i_atom_end = structure.atoms.end();
      t_Atoms :: const_iterator i_atom0 = structure0.atoms.begin();
      for(; i_atom != i_atom_end; ++i_atom, ++i_atom0 )
      {
        atat::rVector3d pos;
        if ( not (i_atom->freeze & t_Atom::FREEZE_X ) )
          { pos[0] = types::t_real(2.0) * (*_i_x ); ++_i_x; }
        else pos[0] = i_atom0->pos[0];
        if ( not (i_atom->freeze & t_Atom::FREEZE_Y ) )
          { pos[1] = types::t_real(2.0) * (*_i_x ); ++_i_x; }
        else pos[1] = i_atom0->pos[1];
        if ( not (i_atom->freeze & t_Atom::FREEZE_Z ) )
          { pos[2] = types::t_real(2.0) * (*_i_x ); ++_i_x; }
        else pos[2] = i_atom0->pos[2];

        i_atom->pos = _strain * pos;
      }

      // Correct stray movements of the center of mass
      __ASSERT(    fixed_index[0] < 0
                or fixed_index[1] < 0
                or fixed_index[2] < 0,
                "fixed_index contains negative indices. Was init() called?\n" )
      __ASSERT(    fixed_index[0] >= structure0.atoms.size()
                or fixed_index[1] >= structure0.atoms.size()
                or fixed_index[2] >= structure0.atoms.size(),
                "fixed_index contains out-of-range indices.\n" )

      types::t_real x =   structure0.atoms[fixed_index[0]].pos[0] 
                        - structure.atoms [fixed_index[0]].pos[0];
      types::t_real y =   structure0.atoms[fixed_index[1]].pos[1] 
                        - structure.atoms [fixed_index[1]].pos[1];
      types::t_real z =   structure0.atoms[fixed_index[2]].pos[2] 
                        - structure.atoms [fixed_index[2]].pos[2];

      if ( Fuzzy::eq(x, 0e0) and Fuzzy::eq(y, 0e0) and Fuzzy::eq(z, 0e0) ) return;
      for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
        { i_atom->pos[0] += x; i_atom->pos[1] += y; i_atom->pos[2] += z; }
    }


    // initializes stuff before minimization
    bool Functional :: init( t_Arg& _arg ) 
    {
      // sets up structure0, needed for fractional vs cartesian shit
      structure0 = structure;

      // Now counts the degrees of freedom
      types::t_unsigned dof = 0;
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_XX ) ) ++dof;
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_XY ) ) ++dof;
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_XZ ) ) ++dof;
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_YY ) ) ++dof;
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_YZ ) ) ++dof;
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_ZZ ) ) ++dof;
      dof += posdofs();

      if( not dof ) return false;
     
      _arg.resize( dof );

      atat::rMatrix3d strain;
      strain.zero(); 
      strain(0,0) = types::t_real(1.0);
      strain(1,1) = types::t_real(1.0);
      strain(2,2) = types::t_real(1.0);

      pack_variables( _arg, strain);
      
      return true;
    }

    types::t_unsigned Functional :: posdofs() 
    {
      fixed_index[0] = -1; fixed_index[1] = -1; fixed_index[2] = -1; 
      types::t_unsigned dof = 0;
      t_Atoms :: iterator i_atom =  structure0.atoms.begin();
      t_Atoms :: iterator i_atom_end =  structure0.atoms.end();
      for( types::t_unsigned n = 0; i_atom != i_atom_end; ++i_atom, ++n ) 
      {
        if ( not (i_atom->freeze & t_Atom::FREEZE_X ) ) ++dof;
        else if (fixed_index[0] == -1 ) fixed_index[0] = n;
        if ( not (i_atom->freeze & t_Atom::FREEZE_Y ) ) ++dof;
        else if (fixed_index[1] == -1 ) fixed_index[1] = n;
        if ( not (i_atom->freeze & t_Atom::FREEZE_Z ) ) ++dof;
        else if (fixed_index[2] == -1 ) fixed_index[2] = n;
      }
      
      if( fixed_index[0] == -1 ) fixed_index[0] = 0;
      if( fixed_index[1] == -1 ) fixed_index[1] = 0;
      if( fixed_index[2] == -1 ) fixed_index[2] = 0;

      return dof;
    }

    // variables is expected to be of sufficient size!!
    // call init() first
    void Functional :: pack_variables( t_Arg& _arg, const atat::rMatrix3d& _strain ) const
    {
      // finally, packs vff format into function::Base format
      t_Arg :: iterator i_var = _arg.begin();
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_XX ) )
        { *i_var = _strain(0,0); ++i_var; }
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_YY ) )
        { *i_var = _strain(1,1); ++i_var; }
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_ZZ ) )
        { *i_var = _strain(2,2); ++i_var; }
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_XY ) )
        { *i_var = 0.5*(_strain(1,0) + _strain(0,1)); ++i_var; }
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_XZ ) )
        { *i_var = 0.5*(_strain(2,0) + _strain(0,2)); ++i_var; }
      if ( not (structure0.freeze & Crystal::Structure::FREEZE_YZ ) )
        { *i_var = 0.5*(_strain(2,1) + _strain(1,2)); ++i_var; }

      pack_positions( i_var );
    }

    void Functional :: pack_positions( t_Arg :: iterator &_i_var ) const
    {
       t_Atoms :: const_iterator i_atom =  structure0.atoms.begin();
       t_Atoms :: const_iterator i_atom_end =  structure0.atoms.end();
       for(; i_atom != i_atom_end; ++i_atom )
       {
         if ( not (i_atom->freeze & t_Atom::FREEZE_X ) )
           { *_i_var = i_atom->pos[0] * types::t_real(0.5); ++_i_var; }
         if ( not (i_atom->freeze & t_Atom::FREEZE_Y ) )
           { *_i_var = i_atom->pos[1] * types::t_real(0.5); ++_i_var; }
         if ( not (i_atom->freeze & t_Atom::FREEZE_Z ) )
           { *_i_var = i_atom->pos[2] * types::t_real(0.5); ++_i_var; }
       }
    }

    void Functional :: pack_gradients( const atat::rMatrix3d& _stress, 
                                       t_GradientArg _grad) const
    {
      t_GradientArg i_grad(_grad);

      // first, external stuff
      if ( not (structure.freeze & Crystal::Structure::FREEZE_XX) )
        *i_grad = _stress(0,0), ++i_grad;
      if ( not (structure.freeze & Crystal::Structure::FREEZE_YY) ) 
        *i_grad = _stress(1,1), ++i_grad;
      if ( not (structure.freeze & Crystal::Structure::FREEZE_ZZ) ) 
        *i_grad = _stress(2,2), ++i_grad;
      if ( not (structure.freeze & Crystal::Structure::FREEZE_XY) ) 
        *i_grad = 0.5 * (_stress(0,1) + _stress(1,0)), ++i_grad;
      if ( not (structure.freeze & Crystal::Structure::FREEZE_XZ) ) 
        *i_grad = 0.5 * (_stress(0,2) + _stress(2,0)), ++i_grad;
      if ( not (structure.freeze & Crystal::Structure::FREEZE_YZ) ) 
        *i_grad = 0.5 * (_stress(1,2) + _stress(2,1)), ++i_grad;

      // then atomic position stuff
      t_Centers :: const_iterator i_center = centers.begin();
      t_Centers :: const_iterator i_end = centers.end();
      t_Atoms :: const_iterator i_atom0 = structure0.atoms.begin();
      i_center = centers.begin();
      for (; i_center != i_end; ++i_center, ++i_atom0)
      {
        const atat::rVector3d& gradient = i_center->gradient;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_X) ) 
          *i_grad = gradient[0], ++i_grad;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_Y) ) 
          *i_grad = gradient[1], ++i_grad;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_Z) ) 
          *i_grad = gradient[2], ++i_grad;
      }
      __MPICODE
      (
        typedef boost::remove_pointer<t_GradientArg>::type t_Type;
        // boost does not do in-place reduction.
        BOOST_MPI_CHECK_RESULT
        (
          MPI_Allreduce,
          (
            MPI_IN_PLACE, _grad, size_t(i_grad - _grad),
            boost::mpi::get_mpi_datatype<t_Type>(*_grad),
            boost::mpi::is_mpi_op< std::plus<t_Type>, t_Type>::op(), (MPI_Comm) MPI_COMM
          )
        );
//       MPI_Comm __commC = (MPI_Comm) ( MPI_COMM ) ;
//       boost::mpi::all_reduce( MPI_COMM, _grad, size_t(i_grad - _grad),
//                               std::plus<types::t_real>());
      )
    }

    void Functional :: gradient( const t_Arg& _arg, t_GradientArg _i_grad ) const
    {
      atat::rMatrix3d strain; strain.zero();
      t_Return energy(0);
      foreach( const t_Center& center, centers ) center.gradient = atat::rVector3d(0,0,0);

      // unpacks variables into vff atomic_center and strain format
      unpack_variables( _arg, strain);

      // computes K0
      atat::rMatrix3d K0 = (!(~strain));

      // computes energy and gradient
      stress.zero();
      LADA_MPI_SPLIT_LOOP( t_Centers :: const_iterator, center, centers, MPI_COMM )
      for (; i_center != i_center_end; ++i_center)
        energy += functionals[i_center->kind()].
                       evaluate_with_gradient( *i_center, strain, stress, K0 );
#     ifdef _MPI
        energy = boost::mpi::all_reduce( MPI_COMM, energy, std::plus<types::t_real>() ); 
        // boost does not do in-place reduction.
        BOOST_MPI_CHECK_RESULT
        (
          MPI_Allreduce,
          (
            MPI_IN_PLACE, &(stress.x[0][0]), 9u,
            boost::mpi::get_mpi_datatype<types::t_real>(stress(0,0)),
            boost::mpi::is_mpi_op< std::plus<types::t_real>, types::t_real>::op(),
            (MPI_Comm) MPI_COMM
          )
        );
#     endif

      // now repacks into function::Base format
      pack_gradients(stress, _i_grad);
    }

  } // namespace vff
} // namespace LaDa
