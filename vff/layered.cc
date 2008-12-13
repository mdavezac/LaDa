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
      return;

      // First, lets create an orthonormal vector to u
      atat::rVector3d a1 = Fuzzy::eq( u(0), 0e0 ) ? 
                             ( Fuzzy::eq( u(1), 0e0 ) ? 
                                atat::rVector3d(1, 0, 0):
                                atat::rVector3d(0, u(2), -u(1) )  ): 
                             atat::rVector3d( -u(2) -u(1), u(0), u(0) );
      a = ( 1.0 / std::sqrt( atat::norm2(a1) ) );
      a1 =  a * a1;

      // Then, lets create another... 
      atat::rVector3d a2;
      a2( 0 ) = u(1) * a1(2) - u(2) * a1(1);
      a2( 1 ) = u(2) * a1(0) - u(0) * a1(2);
      a2( 2 ) = u(0) * a1(1) - u(1) * a1(0);

      // Finally, we create a transition matrix from a "u" basis to a cartesian basis.
      atat::rMatrix3d T; T.set_column(0, u); T.set_column(1, a1); T.set_column(2,a2);
      atat::rMatrix3d S; S.zero(); S(0,0) = 1;
      template_strain = T * S * (~T);
    }

    types::t_real Layered :: evaluate()
    {
      unpack_variables(strain);
      return t_Base::energy();
    }


    // initializes stuff before minimization
    bool Layered :: init()
    {
      if( not is_fixed_by_input ) create_template_strain();
      // sets up structure0, needed for fractional vs cartesian shit
      structure0 = structure;

      // Now counts the degrees of freedom
      types::t_unsigned dof = 1 + posdofs();

      function::Base<> :: resize( dof );
      if ( not variables ) return false;

      strain.zero(); 
      strain(0,0) = 1.0;
      strain(1,1) = 1.0;
      strain(2,2) = 1.0;
      pack_variables(strain);
      
      return true;
    }

    // variables is expected to be of sufficient size!!
    void Layered :: pack_variables( const atat::rMatrix3d& _strain)
    {
      __ASSERT( variables->size() == 0, "Too few variables\n" )
      // finally, packs vff format into function::Base format
      iterator i_var = variables->begin();
      *i_var = u * (_strain * u) - 1.0;
      ++i_var;
      pack_positions( i_var );
    }

    // Unpacks opt::Function_Base::variables into Vff::Layered format
    void Layered :: unpack_variables(atat::rMatrix3d& strain)
    {
      __ASSERT( variables->size() < 3, "Too few variables.\n" )
      std::vector<types::t_real> :: const_iterator i_x = variables->begin();

      strain = (*i_x) * template_strain; ++i_x;
      strain(0,0) += 1.0;
      strain(1,1) += 1.0;
      strain(2,2) += 1.0;

      // compute resulting cell vectors
      structure.cell = strain * structure0.cell;
      unpack_positions( strain, i_x );
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
  } // namespace vff
} // namespace LaDa
