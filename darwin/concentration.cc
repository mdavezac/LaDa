//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "concentration.h"
#include <crystal/fourier.h>

#include <opt/debug.h>       // std::runtime_error

namespace LaDa
{
  types::t_real set_concentration( Crystal::Structure &_str,
                                   types::t_real _target )
  {
    types::t_unsigned N = (types::t_int) _str.atoms.size();
    types::t_complex  *hold = new types::t_complex[ N ];

    __DOASSERT( not hold,
                " Could not allocate memory in set_concentration.\n" )

    types::t_complex  *i_hold = hold;
    Crystal::Fourier( _str.atoms.begin(), _str.atoms.end(),
                       _str.k_vecs.begin(), _str.k_vecs.end(),
                       i_hold );

    Crystal::Structure::t_Atoms::iterator i_atom = _str.atoms.begin();
    Crystal::Structure::t_Atoms::iterator i_atom_end = _str.atoms.end();
    types :: t_int result = 0; 
    i_hold = hold;
    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
    {
      if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T )
        ( i_atom->type > 0 ) ? ++result : --result;
      else 
      {
        ( std::real(*i_hold) > 0 ) ? ++result : --result;
        i_atom->type = ( std::real(*i_hold) > 0 ) ? 1.0: -1.0;
      }
    }

    if ( _target == -2.0 )
    {
      delete[] hold;
      return (double) result  / (double) N;
    }

    types::t_int to_change = static_cast<types::t_int>( (double) N * _target ) - result;
    if (  not to_change ) 
    {
      delete[] hold;
      return _target;
    }

    i_atom = _str.atoms.begin();
    i_hold = hold;
    for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
      if ( i_atom->freeze & Crystal::Structure::t_Atom::FREEZE_T )
        *i_hold = to_change;
      else if ( to_change > 0 and std::real( *i_hold ) > 0 )
        *i_hold = to_change;
      else if ( to_change < 0 and std::real( *i_hold ) < 0 )
        *i_hold = to_change;

    i_atom = _str.atoms.begin();
    types::t_int which;
    do
    {
      i_hold = hold;
      which = 0;
      types::t_real maxr = 1.0;
      if ( to_change > 0 )
      {
        for( types::t_unsigned i = 0; i < N; ++i, ++i_hold ) 
          if (     ( maxr == 1.0 and std::real( *i_hold ) < 0 )
               or  (  std::real( *i_hold ) < 0 and maxr <  std::real( *i_hold ) ) )
          {
            maxr = std::real( *i_hold );
            which = i;
          }
        ( i_atom + which )->type = 1.0;
        *( hold + which ) = to_change;
        result+=2; to_change-=2;
        continue;
      }
      maxr = -1.0;
      for( types::t_unsigned i = 0; i < N; ++i, ++i_hold ) 
        if (     ( maxr == -1.0 and std::real( *i_hold ) > 0 )
             or  (  std::real( *i_hold ) > 0 and maxr >  std::real( *i_hold ) ) )
        {
          maxr = std::real( *i_hold );
          which = i;
        }
      ( i_atom + which )->type = -1.0;
      *( hold + which ) = to_change;
      result-=2; to_change+=2;

    } while (to_change); 
    
    delete[] hold;
    return (double) result  / (double) N;
  }


  bool X_vs_y :: Load ( const TiXmlElement &_node )
  {
    const TiXmlElement *parent;
    std::string str;

    // This whole section tries to find a <Functional type="vff"> tag
    // in _node or its child
    str = _node.Value();
    if ( str.compare("Functional" ) != 0 )
         parent = _node.FirstChildElement("Functional");
    else parent = &_node;

    while (parent)
    {
      str = "";
      if ( parent->Attribute( "type" )  )
        str = parent->Attribute("type");
      if ( str.compare("Concentration" ) == 0 )
        break;
      parent = parent->NextSiblingElement("Functional");
    }
    if ( not parent )
    {
      std::cerr << "Could not find an <Functional type=\"Concentration\"> tag in input file" 
                << std::endl;
      return false;
    } 

    if (     parent->Attribute("x0") 
         and parent->Attribute("y0") )
    {
      double d;
      parent->Attribute("x0", &d); x0 = (types::t_real) d;
      parent->Attribute("y0", &d); y0 = (types::t_real) d;
      if ( x0 > 0 and x0 < 1 and y0 > 0 and y0 < 1 )
      {
        x0 = 2.0*x0 - 1.0;
        y0 = 2.0*y0 - 1.0;
        single_c = true;
        return true;
      }
      std::cerr << "Incorrect values for concentrations " << std::endl
                << " x0 = " << x0 << " y0 = " << y0 << std::endl;
    }
    if (    ( not parent->Attribute("a") )
         or ( not parent->Attribute("b") )
         or ( not parent->Attribute("c") ) )
    {
      std::cerr << "Concentration functional is missing some input "<<std::cerr;
      return false;
    }
    parent->Attribute("a", &a);
    parent->Attribute("b", &b);
    parent->Attribute("c", &c);

    if (    ( b*b - 4*a*c < 0.0 )
         or ( b*b - 4*(a-1.0)*c < 0.0 ) )
    {
      std::cerr << "Equation incorrect on input... should be x = c +b*y + a*y*y "<<std::cerr;
      std::cerr << " with x,y in [0,1] "<<std::cerr;
    }

    // changing to x,y in [-1,1]
    c = 2.0*c + b + a*0.5 - 1.0;
    b += a;
    a *= 0.5;

    if ( parent->Attribute("x0") )
    {
      parent->Attribute("x0", &x0);
      if ( x0 > 0 and x0 < 1 and can_inverse(x0) ) 
      {
        x0 = 2.0*x0 - 1.0;
        y0 = get_y(x0);
        single_c = true;
        return true;
      }
      std::cerr << "Incorrect values for concentrations " << std::endl
                << " x = " << x0 << std::endl;
    }
    if ( parent->Attribute("y0") )
    {
      parent->Attribute("y0", &y0);
      if ( y0 > 0 and y0 < 1 )
      {
        y0 = 2.0*y0 - 1.0;
        x0 = get_x(y0);
        single_c = true;
        return true;
      }
      std::cerr << "Incorrect values for concentrations " << std::endl
                << " y = " << y0 << std::endl;
    }

    if (     std::abs( a ) < types::tolerance 
         and std::abs( b ) < types::tolerance )
    {
      std::cerr << " Equation will run into numerical errors..."
                << " b and a coefficients too small " << std::endl;
      return false;
    }
    return true;
  }

  // verifies that _y leads to x in [-1;1] 
  // expects polynomial to be monotonous between y(1) and y(-1)
  bool X_vs_y :: can_inverse( types::t_real _x)
  {
    types::t_real c0 = c + b + a; // y = 1
    types::t_real c1 = c - b + a; // y = -1
    if ( c0 > c1 )
      return _x <= c0 and _x >= c1;
    return _x >= c0 and _x <= c1;
  }


  std::string X_vs_y :: print () const 
  {
    std::ostringstream sstr;
    if ( single_c )
      sstr << "Single Concentration, x0 = " << x0
           << ", y0 = " << y0;
    else
      sstr << "Concentration function x = " 
           << a << " * y^2 + "  << b << " * y + "  << c;
    return sstr.str();
  }
} // namespace LaDa
