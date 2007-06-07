#ifndef _FITNESS_H_
#define _FITNESS_H_


#include<iostream>
#include "opt_functors.h"

namespace function {

  template<typename T_QUANTITY, typename T_BASELINE>
  class Fitness: public Base<typename T_QUANTITY::t_Type, typename T_QUANTITY::t_Container>
  {
    public:
      typedef T_QUANTITY t_Quantity;
      typedef T_BASELINE t_Baseline;
      typedef typename t_Quantity::t_Type t_Type;
      typedef typename t_Quantity::t_Container t_Container;
      using Base<t_Type, t_Container> :: size;
      using Base<t_Type, t_Container> :: get_concentration;

    protected:
      using Base<t_Type, t_Container> :: variables;

    private: 
      t_Quantity *quantity;
      t_Baseline *baseline;

    public:
      // constructors
      Fitness() 
        { quantity = NULL; baseline = NULL; }
      explicit
        Fitness(t_Quantity * const q, t_Baseline * const b) 
        {
          variables = q->get_variables();
          quantity = q; 
          baseline = b; 
          check_pointers("Fitness(t_Quantity *q, t_Baseline *b)");
        }

      // get, set
      virtual void set_variables ( t_Container* const _var )
        { variables = _var; quantity->set_variables(_var); }
      virtual void destroy_variables()
        {
          quantity->destroy_variables();
          variables = NULL;
        }
      void set_quantity ( t_Quantity * const q )
      {
        quantity = q;
        variables = q->get_variables();
      }
      void set_baseline ( t_Baseline * const b )
        { baseline = b; }
      t_Quantity* get_quantity () const
        { return quantity; }
      t_Baseline* get_baseline () const
        { return baseline; }
    
      // required minimizer behaviors
      t_Type evaluate() 
      { 
        check_pointers("evaluate()");
        return (   quantity->evaluate() 
                 - baseline->evaluate( get_concentration() ) );
      }
      virtual void evaluate_gradient(t_Type* const _i_grad) 
      {
        check_pointers("evaluate_gradient(t_Type *gradient)");
        quantity->evaluate_gradient(_i_grad);
        t_Type N = t_Type(variables->size());
        t_Type x = get_concentration();
        t_Type g_shift = baseline->evaluate_gradient( x ) / N;
        t_Type* i_first = _i_grad, *i_last = _i_grad + size(); 
        for ( ; i_first != i_last; ++i_first )
          *i_first -= g_shift;
      };
      virtual t_Type evaluate_one_gradient(types::t_unsigned _pos)
      {
        t_Type x = get_concentration();
        t_Type value  = quantity->evaluate_one_gradient(_pos);
        t_Type N = t_Type(variables->size());
        t_Type g_shift = baseline->evaluate_gradient( x ) / N;
        return value - g_shift;
      }
      virtual t_Type evaluate_with_gradient(t_Type* const _i_grad)
      {
        t_Type x = get_concentration();
        check_pointers("evaluate_polynome_and_gradient(t_Type *gradient)");
        t_Type N = t_Type(variables->size());
        t_Type value  = quantity->evaluate_with_gradient(_i_grad);
        value -= baseline->evaluate( x );
        t_Type g_shift = baseline->evaluate_gradient( x ) / N;
        t_Type *i_first = _i_grad, *i_last = _i_grad + size(); 
        for ( ; i_first != i_last; ++i_first )
          *i_first -= g_shift;
        return value;
      };

      virtual bool init() { return quantity->init(); }
      virtual bool is_taboo() const
        { return quantity->is_taboo(); }
      // debug
    protected:
      void check_pointers(const char *stream) const
      {
        #ifdef _DEBUG_LADA_
          if ( !quantity or !baseline )
          {
            std::cerr << "Pointers not set in Fitness::" << stream
                      << std::endl;
            exit(0);
          } 
        #endif // _DEBUG_LADA_
      }
  };


} // namespace opt

#endif // _FITNESS_H_
