//
//  Version: $Id$
//
#ifndef _OPT_MINIMIZE_OPT_Ssquared_H_
#define _OPT_MINIMIZE_OPT_Ssquared_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "opt/opt_minimize_opt++.h" // always include before other opt++

// implements S^2 - 1 = 0 type constraints 
// works with package opt++, see include below
#include <NonLinearEquation.h>
#include <NonLinearConstraint.h>
#include <CompoundConstraint.h>
#include <OptQNIPS.h>

#include "analysis/analyze_code.h"

namespace opt {


  template<typename Object> 
  class Minimize_Ssquared : public Minimize_Opt< Object >
  {
    protected:
      using Minimize_Opt<Object>::current_func;
      using Minimize_Opt<Object>::opt_init;
      using Minimize_Opt<Object>::opt_evaluate;
    public:
      inline Minimize_Ssquared() : Minimize_Opt<Object>() {};
      explicit
        inline Minimize_Ssquared(Object* _func) : Minimize_Opt<Object>( _func ) {};
      virtual ~Minimize_Ssquared(){};

      virtual bool minimize()
      {
        START_ANALYSIS( "Minimize_Ssquared :: minimize" )
        types::t_int pb_size = current_func->size();
  
        ColumnVector  x(pb_size), b(pb_size);
        opt_init(pb_size, x);
      
        // construct non-linear constraint object
        NLP constraint_func ( new NLF1 (pb_size, pb_size, opt_evaluate_constraints, 
                              opt_init) );
        b=0.0;
        NonLinearConstraint *s_eqn = new NonLinearEquation( &constraint_func, b, pb_size );
        CompoundConstraint constraint( s_eqn );
        
        //  Create the nonlinear function to minimize
        NLF1 func_to_minimize(pb_size, opt_evaluate, opt_init, &constraint);
        
        // set starting values
        func_to_minimize.setX(x);
        func_to_minimize.eval();
  
        // creates method object
        OptQNIPS method( &func_to_minimize );
        method.setFcnTol(0.1);
        method.setMaxStep(0.01);
        method.setMaxIter(150);
        method.setOutputFile("example2.out", 0);
  
        // optimize!!
        method.optimize();
      
        // set optimal solution to variables
        x = func_to_minimize.getXc();
        func_to_minimize.setX(x);
        func_to_minimize.eval();
        typename Object :: iterator i_var = current_func->begin();
        typename Object :: iterator i_var_end = current_func->end();
        for (types::t_unsigned i = 1; i_var != i_var_end; ++i_var, ++i )
          *i_var = x(i);
        
        method.cleanup();
  
        END_ANALYSIS;
        return true;
      
      } // minimize
      // constraints have all the same form
      // We do not need to specify any structure dependent parameters to evaluate them
      // hence everything is done in the following static function
      // all we do is set current_func->variables = x for consistency
      static void opt_evaluate_constraints(types::t_int mode, const types::t_int ndim, const ColumnVector& x, 
                                           ColumnVector& fx, Matrix &gx, types::t_int &result)
      {
        typename Object :: iterator i_var = current_func->begin();
        typename Object :: iterator i_var_end = current_func->end();
        for (types::t_unsigned i = 1; i_var != i_var_end; ++i_var, ++i )
          *i_var = x(i);
          
        // sets gradients to zero 
        gx = 0.0;
      
        // Each constrain is of the form S_i^2-1
        // the derivatives are 2S_i
        Real *ptr_x = x.Store();
        Real *ptr_fx = fx.Store();
        for (types::t_int i=1; i <= ndim; ++i, ++ptr_x, ++ptr_fx)
        {
          if (mode & NLPGradient)
          {
            gx(i,i) = 2.0 * x(i);
            result = NLPGradient;
          }
          if (mode & NLPFunction)
          {
            fx(i) = x(i)*x(i) - 1.0;
            result = NLPFunction;
          }
        } // for (types::t_int i=1; i < ndim; i++, ptr_x++, ptr_fx++)
      
      } // opt_evaluate_constraints
  };
}


#endif
