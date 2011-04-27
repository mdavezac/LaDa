#include "LaDaConfig.h"

#include <cstdlib>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#ifdef LADA_MPI
# include <boost/mpi/collectives/all_reduce.hpp>
#endif

#include <math/fuzzy.h>
#include <opt/debug.h>
#include <opt/path.h>

#include "layered.h"

namespace LaDa
{
  namespace vff
  { 

    void Layered::create_template_strain()
    {
      // The first vector of the cell should indicate the direction of the
      // layering.
      u = is_fixed_by_input ? direction: structure.cell.col(0);
      template_strain = math::rMatrix3d::Zero(); 
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          template_strain(i,j) = u(i)*u(j);
      const types::t_real b = types::t_real(1.0) / u.squaredNorm();
      const types::t_real a = std::sqrt( b );
      u = a * u;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          template_strain(i,j) *= b;
    }

    Layered::t_Return Layered::operator()(const t_Arg& _arg) const
    {
      math::rMatrix3d strain;
      unpack_variables( _arg, strain );
      return Vff::energy();
    }


    // initializes stuff before minimization
    bool Layered :: init( t_Arg& _arg)
    {
      if( not is_fixed_by_input ) create_template_strain();
      // sets up structure0, needed for fractional vs cartesian shit
      structure0 = structure;

      // Now counts the degrees of freedom
      types::t_unsigned dof = 1 + posdofs();

      _arg.resize( dof );

      math::rMatrix3d strain = math::rMatrix3d::Zero();
      strain(0,0) = types::t_real(1.0);
      strain(1,1) = types::t_real(1.0);
      strain(2,2) = types::t_real(1.0);
      pack_variables( _arg, strain );
      
      return true;
    }

    // variables is expected to be of sufficient size!!
    void Layered :: pack_variables( t_Arg& _arg, const math::rMatrix3d& _strain) const
    {
      // finally, packs vff format into function::Base format
      t_Arg :: iterator i_var = _arg.begin();
      math::rMatrix3d strain = _strain;
      strain(0,0) -= types::t_real(1.0);
      strain(1,1) -= types::t_real(1.0);
      strain(2,2) -= types::t_real(1.0);
      *i_var = 0e0;
      size_t n(0);
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
        {
          if( math::is_zero( template_strain(i,j) ) ) continue;
          *i_var += template_strain(i,j) * strain(i,j);
          ++n;
        }
      if( n ) *i_var /= types::t_real(n);
      ++i_var;
      pack_positions( i_var );
    }

    // Unpacks opt::Function_Base::variables into Vff::Layered format
    void Layered :: unpack_variables( const t_Arg& _arg, math::rMatrix3d& strain ) const
    {
      t_Arg :: const_iterator i_x = _arg.begin();

      strain = (*i_x) * template_strain; ++i_x;
      strain(0,0) += types::t_real(1.0);
      strain(1,1) += types::t_real(1.0);
      strain(2,2) += types::t_real(1.0);

      // compute resulting cell vectors
      structure.cell = strain * structure0.cell;
      unpack_positions( i_x, strain );
    }

    void Layered::pack_gradients(const math::rMatrix3d& _stress, t_GradientArg _grad) const
    {
      t_GradientArg i_grad(_grad);

      // first, external stuff
      *i_grad = 0e0;
      for( size_t i(0); i < 3; ++i )
        for( size_t j(0); j < 3; ++j )
          *i_grad += template_strain(i,j) * _stress(i,j);
      ++i_grad;

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
            boost::mpi::is_mpi_op< std::plus<t_Type>, t_Type>::op(), 
            (MPI_Comm) comm
          )
        );
//       boost::mpi::all_reduce
//       (
//         _comm, 
//         _grad, size_t(i_grad - _grad), _grad,
//         std::plus<types::t_real>()
//       );
      )
    }

    Layered::t_Return Layered::gradient( const t_Arg& _arg, t_GradientArg _i_grad) const
    {
      math::rMatrix3d strain = math::rMatrix3d::Zero();
      t_Return energy(0);
      foreach( const t_Center& center, centers_ ) center.gradient = math::rVector3d(0,0,0);

      // unpacks variables into vff atomic_center and strain format
      unpack_variables(_arg, strain);

      // computes K0
      math::rMatrix3d K0 = (!(~strain));

      // computes energy and gradient
      stress = math::rMatrix3d::Zero();
      LADA_MPI_SPLIT_LOOP( t_Centers :: const_iterator, center, centers_, comm )
      for (; i_center != i_center_end; ++i_center)
        energy += evaluate_center_with_gradient( *i_center, strain, stress, K0 );

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
      // now repacks into function::Base format
      pack_gradients(stress, _i_grad);
      return energy;
    }

  } // namespace vff
} // namespace LaDa
