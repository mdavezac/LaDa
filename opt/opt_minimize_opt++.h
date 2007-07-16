#ifndef _OPT_MINIMIZE_OPTpp_H_
#define _OPT_MINIMIZE_OPTpp_H_

// interface class for minimizing with OPT++
// works with any class Object which contains a number of member
// functions. see opt::function in opt_functors.h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "opt/opt_minimize_base.h"
#ifdef Real // atat vs newmat pb
  #undef Real
#endif


#include <newmat.h>
#include <NLF.h>


namespace opt {


  template<typename Object> 
  class Minimize_Opt : public Minimize_Base< Object >
  {
    protected:
      using Minimize_Base<Object>::current_func;
    public:
      inline Minimize_Opt() : Minimize_Base<Object>(){};
      explicit
        inline Minimize_Opt(Object* _func) : Minimize_Base<Object>( _func ) {};
      virtual ~Minimize_Opt(){};

      static void opt_evaluate(types::t_int mode, const types::t_int ndim, const ColumnVector& x,
                               types::t_real &fx, ColumnVector& gx, types::t_int &result)
      {
        typename Object::iterator i_var = current_func->begin();
        typename Object::iterator i_var_end = current_func->end();
        for (types::t_unsigned i=1; i_var != i_var_end; ++i_var, ++i )
          *i_var = x(i);

        if ( mode & NLPFunction )
        {
          fx = current_func->evaluate();
          result = NLPFunction;
        }
        if ( mode & NLPFunction )
        {
          current_func->evaluate_gradient(gx.Store());
          result = NLPGradient;
        }
      } // opt_evaluate

      static inline void opt_init(const types::t_int ndim, ColumnVector& x)
      {
        typename Object::const_iterator i_var = current_func->begin();
        typename Object::const_iterator i_var_end = current_func->end();
        for (types::t_unsigned i=1; i_var != i_var_end; ++i_var, ++i )
          x(i) = *i_var;
      } // opt_init

  };
        
}


#endif
