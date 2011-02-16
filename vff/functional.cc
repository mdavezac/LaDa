#include "LaDaConfig.h"

#include <cstdlib>

#include <algorithm>
#include <functional>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/filesystem/operations.hpp>
#ifdef LADA_MPI
# include <boost/mpi/collectives/all_reduce.hpp>
#endif

#include <physics/physics.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>

#include "functional.h"
  
namespace LaDa
{
  namespace vff
  { 
    Functional::t_Return Functional::operator()(const t_Arg& _arg) const
    {
      math::rMatrix3d strain;
      unpack_variables( _arg, strain );
      return Vff::energy();
    }

    // Unpacks opt::Function_Base::variables into Vff::Functional format
    void Functional :: unpack_variables(const t_Arg& _arg, math::rMatrix3d& _strain) const
    {
      t_Arg :: const_iterator i_x = _arg.begin();

      for(size_t i(0), freeze(1); i < 3; ++i)
        for(size_t j(0); j < 3; ++j, freeze <<= 1)
          _strain(i,j) = structure.freeze & freeze ? (i==j ? 1:0): (*i_x++);

      // compute resulting cell vectors
      structure.cell = _strain * structure0.cell;
      unpack_positions( i_x, _strain );
    }

    void Functional :: unpack_positions( t_Arg :: const_iterator &_i_x, 
                                         math::rMatrix3d& _strain ) const
    {
      // then computes positions
      t_Atoms :: iterator i_atom = structure.atoms.begin();
      t_Atoms :: iterator i_atom_end = structure.atoms.end();
      t_Atoms :: const_iterator i_atom0 = structure0.atoms.begin();
      for(; i_atom != i_atom_end; ++i_atom, ++i_atom0 )
      {
        math::rVector3d pos;
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
      LADA_NASSERT(    fixed_index[0] < 0
                or fixed_index[1] < 0
                or fixed_index[2] < 0,
                "fixed_index contains negative indices. Was init() called?\n" )
      LADA_NASSERT(    fixed_index[0] >= structure0.atoms.size()
                or fixed_index[1] >= structure0.atoms.size()
                or fixed_index[2] >= structure0.atoms.size(),
                "fixed_index contains out-of-range indices.\n" )

      types::t_real x =   structure0.atoms[fixed_index[0]].pos[0] 
                        - structure.atoms [fixed_index[0]].pos[0];
      types::t_real y =   structure0.atoms[fixed_index[1]].pos[1] 
                        - structure.atoms [fixed_index[1]].pos[1];
      types::t_real z =   structure0.atoms[fixed_index[2]].pos[2] 
                        - structure.atoms [fixed_index[2]].pos[2];

      if ( math::eq(x, 0e0) and math::eq(y, 0e0) and math::eq(z, 0e0) ) return;
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
      for(size_t freeze(1); freeze < 512; freeze <<= 1)
        if( not (structure0.freeze & freeze) ) ++ dof;
      dof += posdofs();

      if( not dof ) return false;
     
      _arg.resize( dof );

      math::rMatrix3d strain;
      strain = math::rMatrix3d::Identity();

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
    void Functional :: pack_variables( t_Arg& _arg, const math::rMatrix3d& _strain ) const
    {
      // finally, packs vff format into function::Base format
      t_Arg :: iterator i_var = _arg.begin();
      for(size_t i(0), freeze(1); i < 3; ++i)
        for(size_t j(0); j < 3; ++j, freeze <<= 1)
        {
          if(structure0.freeze & freeze) continue;
          *i_var = _strain(i,j);
          ++i_var;
        }
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

    void Functional :: pack_gradients(const math::rMatrix3d& _stress, t_GradientArg _grad) const 
    {
      t_GradientArg i_grad(_grad);

      // first, external stuff
      for(size_t i(0), freeze(1); i < 3; ++i)
        for(size_t j(0); j < 3; ++j, freeze <<= 1)
        {
          if(structure.freeze & freeze) continue;
          *i_grad = _stress(i,j);
          ++i_grad;
        }

      // then atomic position stuff
      t_Centers :: const_iterator i_center = centers_.begin();
      t_Centers :: const_iterator i_end = centers_.end();
      t_Atoms :: const_iterator i_atom0 = structure0.atoms.begin();
      i_center = centers_.begin();
      for (; i_center != i_end; ++i_center, ++i_atom0)
      {
        const math::rVector3d& gradient = i_center->gradient;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_X) ) 
          *i_grad = gradient[0], ++i_grad;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_Y) ) 
          *i_grad = gradient[1], ++i_grad;
        if ( not (i_atom0->freeze & t_Atom::FREEZE_Z) ) 
          *i_grad = gradient[2], ++i_grad;
      }
      LADA_MPI_CODE
      (
        typedef boost::remove_pointer<t_GradientArg>::type t_Type;
        // boost does not do in-place reduction.
        BOOST_MPI_CHECK_RESULT
        (
          MPI_Allreduce,
          (
            MPI_IN_PLACE, _grad, size_t(i_grad - _grad),
            boost::mpi::get_mpi_datatype<t_Type>(*_grad),
            boost::mpi::is_mpi_op< std::plus<t_Type>, t_Type>::op(), (MPI_Comm) comm
          )
        );
      )
    }

    Functional::t_Return Functional :: gradient(const t_Arg& _arg, t_GradientArg _i_grad) const
    {
      math::rMatrix3d strain = math::rMatrix3d::Zero();
      t_Return energy(0);
      foreach( const t_Center& center, centers_ ) center.gradient = math::rVector3d(0,0,0);

      // unpacks variables into vff atomic_center and strain format
      unpack_variables( _arg, strain);

      // computes K0
      math::rMatrix3d K0 = strain.transpose().inverse();

      // computes energy and gradient
      stress = math::rMatrix3d::Zero();
      LADA_MPI_SPLIT_LOOP( t_Centers :: const_iterator, center, centers_, comm )
      for (; i_center != i_center_end; ++i_center)
        energy += evaluate_center_with_gradient(*i_center, strain, stress, K0);

      // First repacks into function::Base format
      pack_gradients(stress, _i_grad);
      // Then sums actual results.
#     ifdef LADA_MPI
        energy = boost::mpi::all_reduce( comm, energy, std::plus<types::t_real>() ); 
        // boost does not do in-place reduction.
        BOOST_MPI_CHECK_RESULT
        (
          MPI_Allreduce,
          (
            MPI_IN_PLACE, stress.data(), stress.stride()*3,
            boost::mpi::get_mpi_datatype<types::t_real>(stress(0,0)),
            boost::mpi::is_mpi_op< std::plus<types::t_real>, types::t_real>::op(),
            (MPI_Comm) comm
          )
        );
#     endif
      return energy;
    }

  } // namespace vff
} // namespace LaDa
