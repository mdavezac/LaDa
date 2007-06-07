#ifndef _OPT_MINIMIZE_OPT_WANG_H_
#define _OPT_MINIMIZE_OPT_WANG_H_

// implements -1<= S <= 1 type constraints
// works with package opt++, see include below
#include <opt/opt_minimize_opt++.h> // always include before other opt++
#include <NLP.h>
#include <Constraint.h>
#include <BoundConstraint.h>
#include <OptBCQNewton.h>

namespace opt {


  template<typename Object> 
  class Minimize_Wang : public Minimize_Opt< Object >
  {
    protected:
      using Minimize_Opt<Object>::current_func;
      using Minimize_Opt<Object>::opt_init;
      using Minimize_Opt<Object>::opt_evaluate;

    public:
      inline Minimize_Wang() : Minimize_Opt<Object>() {};
      explicit
        inline Minimize_Wang(Object* _func) : Minimize_Opt<Object>( _func ) {};
      virtual ~Minimize_Wang(){};

      virtual bool minimize()
      {
        START_ANALYSIS( "Minimize_Wang :: minimize" )
        types::t_int pb_size = current_func->size();
      
        ColumnVector  x(pb_size);
        ColumnVector upper(pb_size), lower(pb_size);
        upper = 1.0; lower = 1.0;
        opt_init(pb_size, x);
       
        // construct constraints
        Constraint constraints = new BoundConstraint(pb_size, lower, upper);
        CompoundConstraint cc(constraints);
      
        //  Create the nonlinear function to minimize
        NLF1 func_to_minimize(pb_size, opt_evaluate, opt_init, &cc);
        
        // set starting values
        func_to_minimize.setX(x);
        func_to_minimize.eval();
      
        // creates method object
        OptBCQNewton method( &func_to_minimize );
        method.setFcnTol(0.1);
        method.setMaxStep(0.01);
        method.setMaxIter(150);
      
        // optimize!!
        method.optimize();
      
        // set optimal solution to variables
        x = func_to_minimize.getXc();
        func_to_minimize.setX(x);
        func_to_minimize.eval();
        typename Object :: iterator i_var = current_func->begin();
        typename Object :: iterator i_var_end = current_func->end();
        for (types::t_unsigned i = 1; i_var != i_var_end; ++i_var, ++i )
          ( fabs( x(i) - 1.0 ) < fabs( x(i) + 1.0 ) ) ?
            *i_var = 1.0 : *i_var = -1.0; 
        
        method.cleanup();
      
        END_ANALYSIS;
        return true;
      
      } // minimize
  };
        
}


#endif
