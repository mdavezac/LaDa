#ifndef _CONCENETRATION_H_
#define _CONCENETRATION_H_

#include <stdexcept>       // std::runtime_error
#include "lamarck/structure.h"
#include "opt/types.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

// template < class T_STEP >
// class Concentration
// {
//   public:
//     typedef T_STEP t_step;
//
//   protected:
//     Ising_CE::Structure &structure;
//     types::t_complex *hold;
//   public:
//     types::t_int site;
//
//   public:
//     Concentration( Ising_CE::Structure &_s ) : structure(_s), hold(NULL), site(-1)
//       { hold = new types::t_complex[ _s.atoms.size() ]; }
//     Concentration( const Concentration& _c ) : structure(_c.structure), hold(_c.hold) {};
//     ~Concentration() { delete[] hold; }
//
//     types::t_real operator()( t_step _step = t_step(0.0) ) // add to structure kvectors
//     {
//       types::t_complex  *i_hold = hold;
//       if ( _step != t_step(0.0) ) 
//         structure.k_vecs.front().type += _step;
//       Ising_CE::fourrier_to_rspace( structure.atoms.begin(), structure.atoms.end(),
//                                     structure.k_vecs.begin(), structure.k_vecs.end(),
//                                     i_hold );
//       
//       types :: t_int result = 0; 
//       types::t_unsigned N = 0;
//       i_hold = hold;
//       Ising_CE::Structure::t_Atoms::const_iterator i_atom = structure.atoms.begin();
//       Ising_CE::Structure::t_Atoms::const_iterator i_atom_end = structure.atoms.end();
//
//       if ( site == -1 )
//       {
//         for (; i_atom != i_atom_end; ++i_atom, ++i_hold, ++N)
//           if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) 
//             ( i_atom->type > 0 ) ? ++result : --result;
//           else
//             ( std::real(*i_hold) > 0 ) ? ++result : --result;
//         return (double) result  / (double) N;
//       }
//
//       for (; i_atom != i_atom_end; ++i_atom, ++i_hold, ++N)
//         if ( site == i_atom->site )
//         {
//           if ( i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T )
//             ( i_atom->type > 0 ) ? ++result : --result;
//           else 
//             ( std::real(*i_hold) > 0 ) ? ++result : --result;
//         }
//
//       return (double) result  / (double) N;
//     }
//     void operator()( Ising_CE::Structure &_str) 
//     {
//       types::t_complex  *i_hold = hold;
//       Ising_CE::Structure::t_Atoms::iterator i_atom = structure.atoms.begin();
//       Ising_CE::Structure::t_Atoms::iterator i_atom_end = structure.atoms.end();
//       if ( site == -1 )
//       {
//         for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
//           if ( not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
//             i_atom->type = ( std::real(*i_hold) > 0 ) ? 1.0 : -1.0;
//         return;
//       }
//       for (; i_atom != i_atom_end; ++i_atom, ++i_hold)
//         if (     i_atom->site == site 
//              and not (i_atom->freeze & Ising_CE::Structure::t_Atom::FREEZE_T ) )
//           i_atom->type = ( std::real(*i_hold) > 0 ) ? 1.0 : -1.0;
//     }
// };

  class X_vs_y
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<X_vs_y>( X_vs_y & );
#endif
    protected:
      types::t_real a, b, c;
      types::t_real x0;
      types::t_real y0;
      bool singlec;

    public:
      X_vs_y() : a(0), b(0), c(0), x0(0), y0(0), singlec(false) {}
      X_vs_y( const X_vs_y &_c) : a(_c.a), b(_c.b), c(_c.c), 
                                  x0(_c.x0), y0(_c.y0), singlec(_c.singlec) {}
      bool Load( const TiXmlElement &_node );
      types::t_real get_x( types::t_real _y ) 
      {
        if ( singlec ) return x0;
        return c + b * _y + a * _y * _y; 
      }
      void set_xy( types::t_real _x, types::t_real _y )
      {
        x0 = _x; y0 = _y; singlec = true;
      } 
      types::t_real get_y() { return y0; }
      types::t_real get_x() { return x0; }
      types::t_real get_y( types::t_real _x )
      {
        if ( singlec ) return y0;
        if ( std::abs ( a ) < types::tolerance )
          return ( _x - c ) / b;
       
        types::t_real det = b*b - 4.0 * (c-_x) * a; 
        if ( det < 0 )
        {
          std::cerr << "Error when using Concentration::get_y(" << _x<<")" << std::endl
                    << "determinent is negative, " << det << std::endl;
          throw std::runtime_error("");
        }
        det = std::sqrt(det);
        types::t_real u = 1.0 / ( 2.0 * a );  
        types::t_real r0 =  (-b + det ) * u;
        types::t_real r1 =  (-b - det ) * u;
        if ( std::abs(r0 - 1.0 ) < types::tolerance ) r0 = 1.0;
        if ( std::abs(r0 + 1.0 ) < types::tolerance ) r0 = -1.0;
        if ( std::abs(r1 - 1.0 ) < types::tolerance ) r1 = 1.0;
        if ( std::abs(r1 + 1.0 ) < types::tolerance ) r1 = -1.0;
        if (     ( r0 < -1.0 or r0 > 1.0 )
             and ( r1 < -1.0 or r1 > 1.0 ) )
        {
          std::cerr << a + b +c << " " << a - b + c << std::endl;
          std::cerr << "Error when using Concentration::get_y(" << _x<< ")" << std::endl;
          std::cerr << " r0= " << r0  << " and r1= " << r1 << std::endl;
          throw std::runtime_error("");
        }
        if ( r0 < -1.0 or r0 > 1.0 ) 
          return r1; 
        return r0;
      }
      bool can_inverse( types::t_real _x );
      bool is_singlec () const { return singlec; }

  };

types::t_real set_concentration( Ising_CE::Structure &_str,
                                 types::t_real _target = -2.0);


#endif 
