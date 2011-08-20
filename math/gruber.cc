#include "LaDaConfig.h"

#include "misc.h"

#include "gruber.h"
#include <root_exceptions.h>


namespace LaDa
{
  namespace math
  {
    types::t_int const max_no_change = 2;

    //! Implementation of the Gruber algorithm.
    //! Class implementations simplify the number of arguments to pass.
    struct Gruber
    {
      rMatrix3d cell, rinv; 
      size_t itermax;
      types::t_int nochange, iterations;
      types::t_real a, b, c, d, e, f;
      types::t_real last_a, last_b, last_c; 
      Gruber   (types::t_int _itermax)
             : cell(rMatrix3d::Zero()), itermax(_itermax), nochange(0)  {}

      rMatrix3d operator()(rMatrix3d const &_in) 
      {
        iterations = 0;
        cell = _in;
        rMatrix3d metric = (~cell) * cell;
        a = metric(0,0);
        b = metric(1,1);
        c = metric(2,2);
        d = 2e0 * metric(1,2);
        e = 2e0 * metric(0,2);
        f = 2e0 * metric(0,1);
        rinv = rMatrix3d::Identity(); 
        nochange = 0;
        size_t iter(0); 
        last_a = -a;
        last_b = -b;
        last_c = -c;
        while(step());
        return cell * rinv;
      };

      bool def_gt_0()
      {
        int z, p;
        def_test(z, p);
        return p == 3 or (z==0 and p==1);
      }

      void def_test(int &_zero, int &_positive)
      {
        _zero = 0; _positive = 0;
        if(lt(0e0, d)) _positive += 1;
        else if(eq(d, 0e0)) _zero += 1;
        if(lt(0e0, e)) _positive += 1;
        else if(eq(e, 0e0)) _zero += 1;
        if(lt(0e0, f)) _positive += 1;
        else if(eq(f, 0e0)) _zero += 1;
      }

      bool step()
      {
        if(gt(a, b) or (eq(a, b) and gt(std::abs(d), std::abs(e)))) n1_action();
        if(gt(b, c) or (eq(b, c) and gt(std::abs(e), std::abs(f)))) { n2_action(); return true; }
        if(def_gt_0()) n3_action();
        else
        { 
          n4_action();
          if(not is_changing()) return false;
        }
        if(    gt( std::abs(d), b) 
            or (eq(d, b)  and lt(e+e, f)) 
            or (eq(d, -b) and lt(f, 0e0)) ) { n5_action(); return true; } 
        if(    gt( std::abs(e), a) 
            or (eq(e, a)  and lt(d+d, f)) 
            or (eq(e, -a) and lt(f, 0e0)) ) { n6_action(); return true; } 
        if(    gt( std::abs(f), a) 
            or (eq(f, a)  and lt(d+d, e)) 
            or (eq(f, -a) and lt(e, 0e0)) ) { n7_action(); return true; } 
        if(    lt(d+e+f+a+b, 0e0)
            or (eq(f, a)  and lt(d+d, e)) 
            or (eq(d+e+f+a+b, 0e0) and gt(a+a+e+e+f, 0e0))) { n8_action(); return true; }
        return false;
      }

      bool is_changing()
      {
        if(    details::no_opt_change_test(a, last_a)
            || details::no_opt_change_test(b, last_b)
            || details::no_opt_change_test(c, last_c) )
         nochange = 0; 
        else ++nochange;
        last_a = a;
        last_b = b;
        last_c = c;
        return nochange < max_no_change;
      }

      void cb_update( types::t_real a00, types::t_real a01,  types::t_real a02, 
                      types::t_real a10, types::t_real a11,  types::t_real a12, 
                      types::t_real a20, types::t_real a21,  types::t_real a22 )
      {
        if(itermax != 0 and iterations >= itermax) BOOST_THROW_EXCEPTION(error::stop_iteration());
        rMatrix3d update;
        update << a00, a01, a02, a10, a11, a12, a20, a21, a22;
        rinv *= update;
        ++iterations;
      }

      void n1_action()
      {
        cb_update(0, -1, 0, -1, 0, 0, 0, 0, -1);
        std::swap(a, b);
        std::swap(d, e);
      }

      void n2_action()
      {
        cb_update(-1, 0, 0, 0, 0, -1, 0, -1, 0);
        std::swap(b, c);
        std::swap(e, f);
      }

      void n3_action()
      {
        types::t_int i(1), j(1), k(1);
        if (lt(d, 0e0)) i = -1;
        if (lt(e, 0e0)) j = -1;
        if (lt(f, 0e0)) k = -1;
        cb_update(i, 0, 0, 0, j, 0, 0, 0, k);
        d = std::abs(d);
        e = std::abs(e);
        f = std::abs(f);
      }
      void n4_action()
      {
        iVector3d vec(1,1,1);
        size_t z(-1);
        if(gt(d, 0e0)) vec(0) = -1; else if( eq(d, 0e0) ) z = 0; 
        if(gt(e, 0e0)) vec(1) = -1; else if( eq(e, 0e0) ) z = 1; 
        if(gt(f, 0e0)) vec(2) = -1; else if( eq(f, 0e0) ) z = 2; 
        if( vec(0) * vec(1) * vec(2) < 0)
        {
          if(z < 0) BOOST_THROW_EXCEPTION(error::internal());
          vec(z) = -1;
        }
        cb_update(vec(0), 0, 0, 0, vec(1), 0, 0, 0, vec(2));
        d = -std::abs(d);
        e = -std::abs(e);
        f = -std::abs(f);
      }

      void n5_action()
      {
        if(d > 0)
        {
          cb_update(1,0,0,0,1,-1,0,0,1);
          c += b - d;
          d -= b + b;
          e -= f;
        }
        else
        {
          cb_update(1,0,0,0,1,1,0,0,1);
          c += b + d;
          d += b + b;
          e += f;
        }
      }
      void n6_action()
      {
        if(e > 0)
        {
          cb_update(1,0,-1,0,1,0,0,0,1);
          c += a - e;
          d -= f;
          e -= a + a;
        }
        else
        {
          cb_update(1,0,1,0,1,0,0,0,1);
          c += a + e;
          d += f;
          e += a + a;
        }
      }
      void n7_action()
      {
        if(f > 0)
        {
          cb_update(1,-1,0,0,1,0,0,0,1);
          b += a - f;
          d -= e;
          f -= a + a;
        }
        else
        {
          cb_update(1,1,0,0,1,0,0,0,1);
          b += a + f;
          d += e;
          f += a + a;
        }
      }
      void n8_action()
      {
        cb_update(1,0,1,0,1,1,0,0,1);
        c += a + b + d + e + f;
        d += b + b + f;
        e += a + a + f;
      }
    };

    //! Computes the Niggli cell of a lattice.
    rMatrix3d gruber(rMatrix3d const &_in, size_t itermax)
    {
      return Gruber(itermax)(_in); 
      try { return Gruber(itermax)(_in); }
      catch(error::stop_iteration &_e) { return rMatrix3d::Zero(); }
    }

  }
}
