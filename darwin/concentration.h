#ifndef _CONCENETRATION_H_
#define _CONCENETRATION_H_

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

    public:
      X_vs_y() : a(0), b(0), c(0) {}
      X_vs_y( const X_vs_y &_c) : a(_c.a), b(_c.b), c(_c.c) {}
      bool Load( const TiXmlElement &_node );
      bool evaluate( types::t_real _x ) 
        { return a + b * _x + c * _x * _x; }
      bool inverse( types::t_real _y )
      {
        if ( std::abs(a - _y) < types::tolerance  )
          { std::cerr << "Error when using Concentration::inverse" <<std::endl; exit(0); }
        types::t_real det = b*b - 4.0 * (a-_y) * c; 
        if ( det < 0 )
          { std::cerr << "Error when using Concentration::inverse" <<std::endl; exit(0); }
        types::t_real u = 1.0 / ( 2.0 * (a - _y ) );  
        types::t_real x0 =  - (b + det ) * u;
        types::t_real x1 =  - (b - det ) * u;
        if (     ( x0 < -1.0 or x0 > 1.0 )
             and ( x1 < -1.0 or x1 > 1.0 ) )
          { std::cerr << "Error when using Concentration::inverse" <<std::endl; exit(0); }
        if ( x0 < -1.0 or x0 > 1.0 ) 
          return x1; 
        return x0;
      }
  };

types::t_real set_concentration( Ising_CE::Structure &_str,
                                 types::t_real _target = -2.0);


#endif 
