//
//  Version: $Id$
//
#ifndef _OBJECTIVE_IMPL_H_
#define _OBJECTIVE_IMPL_H_

namespace Objective
{
  template< class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
  void Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: init( const t_Individual & _indiv)
    { Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = &_indiv; }
  template<class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
    const typename Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: t_Individual* 
      Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: current_indiv = NULL;
  template<class T_GA_TRAITS, class T_QUANTITY_TRAITS, class T_VA_TRAITS >
    typename Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: t_Fitness 
      Base<T_GA_TRAITS, T_QUANTITY_TRAITS, T_VA_TRAITS> :: fitness;


  template< class T_GA_TRAITS >
  inline void Maximize<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_val,
                                                          t_QuantityGradients &_grad,
                                                          t_VA_Type *_i_grad)
  { 
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;
    for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
      *i_grad_result -= *i_grad;
  }
  template< class T_GA_TRAITS >
  inline typename Maximize<T_GA_TRAITS>>::t_ScalarQuantity
    Maximize<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                     t_QuantityGradients& _grad,
                                                     t_VA_Type *_i_grad)  
  { 
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;
    for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
      *i_grad_result -= *i_grad;
    return -_val;
  }

  template< class T_GA_TRAITS >
  inline void Minimize<T_GA_TRAITS> :: evaluate_gradient( const t_Quantity &_val,
                                                          t_QuantityGradients &_grad,
                                                          t_VA_Type *_i_grad)
  { 
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;
    for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
      *i_grad_result += *i_grad;
  }

  template< class T_GA_TRAITS >
  inline typename Minimize<T_GA_TRAITS>>::t_ScalarQuantity
    Minimize<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                     t_QuantityGradients& _grad,
                                                     t_VA_Type *_i_grad)  
  { 
    typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
    typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
    t_VA_Type *i_grad_result = _i_grad;
    for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
      *i_grad_result += *i_grad;
    return _val;
  }
}
#endif
