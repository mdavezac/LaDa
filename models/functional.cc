
//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <opt/debug.h>
#include <algorithm>
#include <boost/lambda/lambda.hpp>

#include "functional.h"
  
namespace LaDa
{
  namespace Models
  { 
    bool Functional :: init( t_Arg& _arg ) 
    {
      forces = structure;
      Clj::t_Arg::t_Atoms :: iterator i_atom = structure.atoms.begin();
      Clj::t_Arg::t_Atoms :: iterator i_atom_end = structure.atoms.end();
      cell0_ = structure.cell;
      scaling_ = 0e0;
      pack_variables( _arg );
      return true; 
    }

    Functional :: t_Return Functional :: operator()( const t_Arg& _arg ) const
    {
      unpack_variables( _arg );
      return Clj::operator()( structure, forces );
    }

    // Unpacks opt::Function_Base::variables into Vff::Functional format
    void Functional :: unpack_variables(const t_Arg& _arg) const
    {
      t_Arg :: const_iterator i_x = _arg.begin();

      if( relaxation == Relaxation::default_ )
      {
        for( size_t i(0); i < 3; ++i )
          for( size_t j(0); j < 3; ++j, ++i_x )
            structure.cell(i,j) = (*i_x);
      }
      else if( relaxation == Relaxation::volume )
      {
        scaling_ = (*i_x++);
        structure.cell = (1e0 + scaling_) * cell0_;
      }

      // then computes positions
      Clj::t_Arg::t_Atoms :: iterator i_atom = structure.atoms.begin();
      Clj::t_Arg::t_Atoms :: iterator i_atom_end = structure.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom )
        for( size_t i(0); i < 3; ++i, ++i_x )  i_atom->pos[i] = *i_x;
    }

    void Functional :: pack_variables( t_Arg& _arg ) const
    {
      _arg.clear();
      if( relaxation == Relaxation::default_ )
      {
        for( size_t i(0); i < 3; ++i )
          for( size_t j(0); j < 3; ++j )
            _arg.push_back( structure.cell(i,j) );
      }
      else if( relaxation == Relaxation::volume ) _arg.push_back( scaling_ );

      Clj::t_Arg::t_Atoms :: iterator i_atom = structure.atoms.begin();
      Clj::t_Arg::t_Atoms :: iterator i_atom_end = structure.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom )
        for( size_t i(0); i < 3; ++i ) _arg.push_back(i_atom->pos[i]);
    }

    void Functional :: pack_gradients( t_GradientArg _grad) const
    {
      t_GradientArg i_grad(_grad);

      // first, external stuff
      if( relaxation == Relaxation::default_ )
      {
        Eigen::Matrix3d const stress( -forces.cell * (~(!structure.cell)) );
        for( size_t i(0); i < 3; ++i )
          for( size_t j(0); j < 3; ++j, ++i_grad )
            *i_grad = stress(i,j);
      }
      else if( relaxation == Relaxation::volume ) 
      {
        Eigen::Matrix3d const stress( -forces.cell * (~(!structure.cell)) );
        *i_grad = 0e0;
        for( size_t i(0); i < 3; ++i )
          for( size_t j(0); j < 3; ++j )
            *i_grad += cell0_(i,j) * stress(i,j);
        ++i_grad;
      }
      Eigen::Matrix3d const matrix( -structure.scale * (~structure.cell) );
      Clj::t_Arg::t_Atoms :: iterator i_atom = forces.atoms.begin();
      Clj::t_Arg::t_Atoms :: iterator i_atom_end = forces.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom )
      {
        Eigen::Vector3d const pos( matrix * i_atom->pos );
        for( size_t i(0); i < 3; ++i, ++i_grad ) *i_grad = pos[i];
      }
    }

    void Functional :: gradient( const t_Arg& _arg, t_GradientArg _i_grad ) const
    {
      unpack_variables( _arg );
      namespace bl = boost::lambda;

      Clj::operator()( structure, forces );
      
      pack_gradients(_i_grad);
    }

  } // namespace vff
} // namespace LaDa
