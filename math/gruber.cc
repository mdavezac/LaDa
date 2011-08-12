#include "LaDaConfig.h"

#include "misc.h"

#include "gruber.h"


namespace LaDa
{
  namespace math
  {
    types::t_real const multiplier = 16e0;
    types::t_int const change = 2;

    //! Computes the Niggli cell of a lattice.
    rMatrix3d gruber(rMatrix3d const &_in, types::t_int itermax = -1)
    {
    }

    //! Implementation of the Gruber algorithm.
    //! Class implementations simplify the number of arguments to pass.
    struct Gruber
    {
      rMatrix3d cell, rinv; 
      types::t_int itermax, nochange;
      types::t_real a, b, c, d, e, f;
      Gruber   (types::t_int _itermax)
             : cell(_in), itermax(_itermax), nochange(0)  {}

      rMatrix3d operator()(rMatrix3d const &_in) 
      {
        rMatrix3d metric = (~cell) * cell;
        a = metric(0,0);
        b = metric(1,1);
        c = metric(2,2);
        d = 2e0 * metric(1,2);
        e = 2e0 * metric(0,2);
        f = 2e0 * metric(0,1);
        rinv = rMatrix::Identity(); 
        nochange(0);
        size_t iter(0); 
        while(step() and ((itermax > 0 and iter < itermax) or itermax <= 0) ) { ++iter; }
        return cell * rinv;
      };

      void def_test(int &_zero, int &_positive)
      {
        _zero = 0; _positive = 0;
        if(gt(d, 0)) _positive += 1;
        else if(le(d, 0)) _zero += 1;
        if(gt(e, 0)) _positive += 1;
        else if(le(e, 0)) _zero += 1;
        if(gt(f, 0)) _positive += 1;
        else if(le(f, 0)) _zero += 1;
      }

      bool step()
      {
        if(le(b, a)) n1_action();
        if(le(c, b)) { n2_action(); return true; }
        int zero, positive;
        if( n3_action() ) return false;
        if( b2_action() ) return true;
        if( b3_action() ) return true;
        if( b4_action() ) return true;
        if( b5_action() ) return true;
        return false;
      }

      void is_changing(
    };

    
  }
}

#endif

