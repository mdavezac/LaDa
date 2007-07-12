#ifndef _BISECTION_H_
#define _BISECTION_H_

namespace method
{
  template<class OPERATION>
  class Bisection
  {
    public:
      typedef OPERATION t_operation;

    protected:
      t_operation &operation;
      types::t_unsigned max_iter;
      types::t_real tolerance;

    public:
      explicit
        Bisection   ( t_operation &_operation )  
                  : operation(_operation), max_iter(32), tolerance(types::tolerance) {};
      ~Bisection() {}

      types::t_real operator()( types::t_real _step,
                                types::t_real _target = 0.0,
                                types::t_unsigned _maxiter = 0,
                                types::t_real _tol = 0)
      {
        if ( _maxiter == 0 ) _maxiter = max_iter;
        if ( _tol <= 0 ) _tol = tolerance;
        types::t_real value = operation( 0.0 ) - _target;
        if ( std::abs( value ) < _tol ) 
          return value;
        types::t_unsigned n = 0;
        bool direction = (value > 0.0);
        bool old_direction = direction;
        types::t_real interval(_step);

        // first finds interval length
//       std::cout << "find interval length, x= " << value + _target << std::endl;
        while (     direction == old_direction 
                and n < _maxiter )
        {

          if ( direction != old_direction )
            break;
//       std::cout << "x = " << value + _target << "  interval: " << interval << std::endl;

          interval += _step;
          value = operation( direction ? -_step: _step ) - _target;

          direction = (value > 0.0);

        }
        if ( n >= _maxiter )
        {
          std::cerr << "Could not find interval length for method::Bisection"<< std::endl;
          return value;
        }
 //      std::cout << "x = " << value + _target << "  interval: " << interval << std::endl
 //                << "now bisect" << std::endl;
        _step = interval;
        n=0;

        // then bisects intervall
        do
        {
          direction = (value > 0.0);
          _step /= 2.0;

          value = operation( direction ? -_step: _step ) - _target;
//         std::cout << "x = " << value + _target << " " << _step << " " << direction << std::endl;

          ++n;
        }
        while (     n < _maxiter  
                and std::abs( value ) > _tol );

//       std::cout << "End : "<< value + _target << std::endl << std::endl;
        return value;
      }
  };
}


#endif
