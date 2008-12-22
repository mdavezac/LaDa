//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdlib>

#include <opt/fuzzy.h>
#include <opt/debug.h>

#include "layered.h"

namespace LaDa
{
  namespace Vff
  { 

    void Layered::create_template_strain()
    {
      // The first vector of the cell should indicate the direction of the
      // layering.
      u = is_fixed_by_input ? direction: structure.cell.get_column(0);
      types::t_real a = 1.0 / std::sqrt( atat::norm2(u) );
      u = a * u;
      template_strain.zero(); 
      template_strain(0,0) = u(0);
      template_strain(1,1) = u(1);
      template_strain(2,2) = u(2);
      std::cout << "template strain: " << template_strain << "\n"; 
      return;
    }

    Layered :: t_Return Layered :: operator()( const t_Arg& _arg ) const
    {
      atat::rMatrix3d strain;
      unpack_variables( _arg, strain );
      t_Return energy = Vff::energy();
      std::cout << "E: " << energy / 16.0217733 << "\n";
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

      atat::rMatrix3d strain; strain.zero(); 
      strain(0,0) = 1.0;
      strain(1,1) = 1.0;
      strain(2,2) = 1.0;
      pack_variables( _arg, strain );
      
      return true;
    }

    // variables is expected to be of sufficient size!!
    void Layered :: pack_variables( t_Arg& _arg, const atat::rMatrix3d& _strain) const
    {
      // finally, packs vff format into function::Base format
      t_Arg :: iterator i_var = _arg.begin();
      *i_var = u * (_strain * u) - 1.0;
      ++i_var;
      pack_positions( i_var );
    }

    // Unpacks opt::Function_Base::variables into Vff::Layered format
    void Layered :: unpack_variables( const t_Arg& _arg, atat::rMatrix3d& strain ) const
    {
      t_Arg :: const_iterator i_x = _arg.begin();

      std::cout << "strain : " << *i_x << "\n";
      strain = (*i_x) * template_strain; ++i_x;
      strain(0,0) += 1.0;
      strain(1,1) += 1.0;
      strain(2,2) += 1.0;

      std::cout << strain << "\n";
      // compute resulting cell vectors
      structure.cell = strain * structure0.cell;
      unpack_positions( i_x, strain );
    }

    bool Layered :: Load_( const TiXmlElement &_node )
    {
      if( not t_Base :: Load( _node ) ) return false;
      if( not _node.Attribute("direction") )
      {
        is_fixed_by_input = false;
        return true;
      }
      bool couldload = true;
      std::istringstream sstr; sstr.str( _node.Attribute("direction") );
      sstr >> direction[0]; if ( sstr.fail() ) couldload = false;
      sstr >> direction[1]; if ( sstr.fail() ) couldload = false;
      sstr >> direction[2]; if ( sstr.fail() ) couldload = false;

      if ( atat::norm2( direction ) < types::tolerance ) couldload = false;

      if ( not couldload )
      {
        is_fixed_by_input = false;
        std::cerr << "Found direction tag, but could not read it correctly.\n";
        return false; 
      }
       
      is_fixed_by_input = true;
      create_template_strain();
      return true;
    }

    bool Layered :: Load( const TiXmlElement &_node )
    {
      const TiXmlElement* parent = opt::find_functional_node( _node, "vff" );
      if( not t_Base :: Load( *parent ) ) return false;
      return Load_( *parent );
    }


     void Layered :: print_out( std::ostream &stream ) const
     {
       t_Base :: print_out( stream );
       if( is_fixed_by_input )
         stream << "Epitaxial Direction fixed on input: " << direction << "\n";
       else stream << "Epitaxial Direction fixed by unit cell\n";
     }

    void Layered :: pack_gradients(const atat::rMatrix3d& _stress, 
                                   t_GradientArg &_grad) const
    {
      t_GradientArg i_grad(_grad);

      // first, external stuff
      *i_grad = u * ( _stress(0,0) * u );
      ++i_grad;

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
    }

    void Layered :: gradient( const t_Arg& _arg, t_GradientArg &_i_grad ) const
    {
      atat::rMatrix3d strain; strain.zero();
      t_Return energy = 0;
      foreach( const t_Center& center, centers ) center.gradient = atat::rVector3d(0,0,0);

      // unpacks variables into vff atomic_center and strain format
      unpack_variables(_arg, strain);

      // computes K0
      atat::rMatrix3d K0 = (!(~strain));

      // computes energy and gradient
      t_Centers :: const_iterator i_center = centers.begin();
      t_Centers :: const_iterator i_end = centers.end();
      stress.zero();
      for (; i_center != i_end; ++i_center)
        energy += functionals[i_center->kind()].
                        evaluate_with_gradient( *i_center, strain, stress, K0 );
      std::cout << "E: " << energy / 16.0217733 << "\n";

      // now repacks into function::Base format
      pack_gradients(stress, _i_grad);
    }

  } // namespace vff
} // namespace LaDa
