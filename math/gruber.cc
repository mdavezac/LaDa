#include "LaDaConfig.h"

#include "misc.h"

#include "gruber.h"
#include <root_exceptions.h>


namespace LaDa
{
  namespace math
  {
    types::t_real const multiplier = 16e0;
    types::t_int const max_no_change = 2;

    //! Implementation of the Gruber algorithm.
    //! Class implementations simplify the number of arguments to pass.
    struct Gruber
    {
      rMatrix3d cell, rinv; 
      types::t_int itermax, nochange, iterations;
      types::t_real a, b, c, d, e, f;
      types::t_real last_a, last_b, last_c; 
      Gruber   (types::t_int _itermax)
             : cell(rMatrix3d::Zero()), itermax(_itermax), nochange(0)  {}

      rMatrix3d operator()(rMatrix3d const &_in) 
      {
        iterations = 0;
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
        if(geq(d, 0e0)) _positive += 1;
        else if(le(d, 0e0)) _zero += 1;
        if(geq(e, 0e0)) _positive += 1;
        else if(le(e, 0e0)) _zero += 1;
        if(geq(f, 0e0)) _positive += 1;
        else if(le(f, 0e0)) _zero += 1;
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

      bool is_changing()
      {
        if(    details::no_opt_change_test(a, last_a, multiplier)
            || details::no_opt_change_test(b, last_b, multiplier)
            || details::no_opt_change_test(c, last_c, multiplier) )
         nochange = 0; 
        else ++nochange;
        last_a = a;
        last_b = b;
        last_c = c;
        return nochange < max_no_change;
      }

      void cb_update(rMatrix3d &_cell)
      {
        if(iterations != 0 and iterations >= itermax) BOOST_THROW_EXCEPTION(error::stop_iteration());
        rinv *= _cell;
        ++iterations;
      }

      void n1_action()
      {
        rMatrix3d update;
        update << 0, -1, 0, -1, 0, 0, 0, 0, -1;
        cb_update(update);
        std::swap(a, b);
        std::swap(d, e);
      }

      void n2_action()
      {
        rMatrix3d update;
        update << -1, 0, 0, 0, 0, -1, 0, -1, 0;
        cb_update(update);
        std::swap(b, c);
        std::swap(e, f);
      }

      bool n3_action()
      {
        if(def_gt_0())
        {
          types::t_int i(1), j(1), k(1);
          if(d < 0) i = -1;
          if(e < 0) j = -1;
          if(f < 0) k = -1;
          rMatrix3d update;
          update << i, 0, 0, 0, j, 0, 0, 0, k;
          cb_update(update);
          d = std::abs(d);
          e = std::abs(e);
          f = std::abs(f);
          return false;
        }
        else
        {
          iVector3d vec(1,1,1);
          size_t z(-1);
          if(d > 0) { vec(0) = -1; z = 0; }
          if(e > 0) { vec(1) = -1; z = 1; }
          if(f > 0) { vec(2) = -1; z = 2; }
          if( vec(0) * vec(1) * vec(2) < 0)
          {
            if(z < 0) BOOST_THROW_EXCEPTION(error::internal());
            vec(z) = -1;
          }
          rMatrix3d update;
          update << vec(0), 0, 0, 0, vec(1), 0, 0, 0, vec(2);
          cb_update(update);
          d = -std::abs(d);
          e = -std::abs(e);
          f = -std::abs(f);
          return true;
        }
      }

      bool b2_action()
      {
        if( not (b < std::abs(d)) ) return false;
        types::t_int j = entier( (d +b) / (2 * b) );
        if(j == 0) return false;
        rMatrix3d update;
        update << 1,0,0, 0,1,-j, 0,0,1;
        cb_update(update);
        c += j * j * b - j*d;
        d -= 2*j*b;
        e -= j*f;
        if(c > 0) BOOST_THROW_EXCEPTION(error::internal());
        return true;
      }

      bool b3_action()
      {
        if( not (a < std::abs(e)) ) return false;
        types::t_int j = entier( (a +e) / (2 * a) );
        if(j == 0) return false;
        rMatrix3d update;
        update << 1,0,-j, 0,1,0, 0,0,1;
        cb_update(update);
        c += j * j * a - j*e;
        d -= j*f;
        e -= 2*j*a;
        if(c > 0) BOOST_THROW_EXCEPTION(error::internal());
        return true;
      }

      bool b4_action()
      {
        if( not (a < std::abs(f)) ) return false;
        types::t_int j = entier( (a +f) / (2 * a) );
        if(j == 0) return false;
        rMatrix3d update;
        update << 1,-j,0, 0,1,0, 0,0,1;
        cb_update(update);
        b += j * j * a - j*f;
        d -= j*e;
        e -= 2*j*a;
        if(b > 0) BOOST_THROW_EXCEPTION(error::internal());
        return true;
      }

      bool b5_action()
      {
        types::t_real const de = d + e;
        types::t_real const fab = f+a+b;
        if( not (std::abs(de+fab) < 0) ) return false;
        types::t_int j = entier( (de +fab) / (2 * fab) );
        if(j == 0) return false;
        rMatrix3d update;
        update << 1,0,-j, 0,1,-j, 0,0,1;
        cb_update(update);
        c += j * j * fab - j*de;
        d -= j*(2*e+f);
        e -= j*(2*a+f);
        if(c > 0) BOOST_THROW_EXCEPTION(error::internal());
        return true;
      }
    
      types::t_int entier(types::t_real const &_ref)
      {
        types::t_int result = static_cast<types::t_int>(_ref);
        if(_ref - result < 0) --result;
        if(not (_ref - result) < 1) ++result;
        return result;
      }
    };

    //! Computes the Niggli cell of a lattice.
    rMatrix3d gruber(rMatrix3d const &_in, types::t_int itermax = -1)
    {
      try { return Gruber(itermax)(_in); }
      catch(error::stop_iteration &_e) { return rMatrix3d::Zero(); }
    }

  }
}
