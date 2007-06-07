#ifndef _DARWIN_MINIMIZER_GENOP_H_
#define _DARWIN_MINIMIZER_GENOP_H_

namespace darwin 
{

#include <eo/eoGenOp.h>

#include <opt/opt_minimize_base.h>
#include <opt/opt_function_base.h>
#include "results.h"

  template<class T_INDIVIDUAL, class T_EVALUATOR, class T_POPULATION>
  class Minimizer_Functional : public function::Base
                               < typename T_EVALUATOR :: t_Functional :: t_Type,
                                 typename T_EVALUATOR :: t_Functional :: t_Container >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR t_Evaluator;
      typedef T_POPULATION t_Population;
      typedef typename t_Evaluator :: t_Functional t_Functional;
      typedef typename T_EVALUATOR :: t_Functional :: t_Type t_Type;
      typedef typename T_EVALUATOR :: t_Functional :: t_Container t_Container;
      typedef typename t_Individual::t_Object t_Object;
      typedef Results<t_Individual, t_Evaluator, t_Population> t_Results;

    protected:
      using function::Base<t_Type, t_Container> :: variables;

    protected:
      t_Results &results;
     
    public:
      Minimizer_Functional( t_Results &_r ) : results(_r) {};
      t_Type evaluate()
      {
        results.invalidate();
        return results.evaluate();
      }
      t_Type evaluate_with_gradient( t_Type *_i_grad )
        { results.invalidate(); return results.evaluate_with_gradient( _i_grad ); }
      void evaluate_gradient( t_Type *_i_grad )
        { return results.evaluate_gradient( _i_grad ); }
      t_Type evaluate_one_gradient( types::t_unsigned _pos )
        { return results.evaluate_one_gradient( _pos ); }
      bool is_taboo() const
        { return results.is_taboo(); }
      void init( t_Individual &_indiv ) 
      {
        results.init( _indiv ); 
        variables = results.get_objective_variables();
      }
      void set_object( t_Individual &_indiv ) 
        { results.set_object( _indiv ); }
      bool init() { return true; } // initialization is done above
        
  };

  template<class T_INDIVIDUAL, class T_EVALUATOR, class T_POPULATION>
  class MinimizerGenOp : public eoGenOp<T_INDIVIDUAL>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR t_Evaluator;
      typedef T_POPULATION t_Population;
      typedef typename t_Evaluator :: t_Functional t_Functional;
      typedef typename minimizer::Base< t_Functional > t_Minimizer;
      typedef typename t_Individual::t_Object t_Object;
      typedef Results<t_Individual, t_Evaluator, t_Population> t_Results;
      typedef Minimizer_Functional<t_Individual, t_Evaluator, t_Population> t_MFunctional;

    protected:
      t_Minimizer &minimizer;
      t_MFunctional functional;

    public:
      explicit
        MinimizerGenOp   ( t_Minimizer &_m, t_Results &_r )
                       : minimizer(_m), functional( _r ) {};

      unsigned max_production(void) { return 1; } 
   
      void apply(eoPopulator<t_Individual>& _pop)
      {
        functional.init( *_pop );
        minimizer( functional );
        functional.evaluate();
        (*_pop).invalidate();
      }
      virtual std::string className() const {return "Darwin::MinimizerGenOp";}
  };


} // namespace darwin

#endif // _DARWIN_MINMIZER_GENOP_H_
