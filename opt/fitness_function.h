#ifndef _FITNESS_H_
#define _FITNESS_H_

#include "LaDaConfig.h"

#include<iostream>

#include "function_functors.h"

namespace LaDa
{
  namespace function {

    //! \brief Creates a functor which is the distance from a quantity to a baseline.
    //! \details A quantity is any class derived of function::Base. A baseline is
    //!          any class which takes a scalar as an argument. Bothe quantity
    //!          and baseline must return a scalar and implement the usual
    //!          routines for a function::Base derived type. The baseline takes
    //!          as an argument the average of the components of the argument of
    //!          the quantity.
    template<typename T_QUANTITY, typename T_BASELINE>
    class Fitness: public Base<typename T_QUANTITY::t_Type, typename T_QUANTITY::t_Container>
    {
      public:
        //! The type of the quantity
        typedef T_QUANTITY t_Quantity;
        //! The type of the baseline
        typedef T_BASELINE t_Baseline;
        //! see function::Base::t_Type
        typedef typename t_Quantity::t_Type t_Type;
        //! see function::Base::t_Container
        typedef typename t_Quantity::t_Container t_Container;

      public:
        using Base<t_Type, t_Container> :: size;
        using Base<t_Type, t_Container> :: get_concentration;

      protected:
        using Base<t_Type, t_Container> :: variables;

      private: 
        //! pointer to the quantity.
        t_Quantity *quantity;
        //! pointer to the baseline.
        t_Baseline *baseline;

      public:
        //! Constructor.
        Fitness() : quantity(NULL), baseline(NULL) {}
        //! Constructor and Initializer.
        explicit Fitness(t_Quantity * const q, t_Baseline * const b);

        //! Sets the variables of the quantity
        virtual void set_variables ( t_Container* const _var )
          { variables = _var; quantity->set_variables(_var); }
        //! Destroy the variables of the quantity
        virtual void destroy_variables();
        //! Sets the pointer to the quantity
        void set_quantity ( t_Quantity * const q );
        //! Sets the pointer to the baseline
        void set_baseline ( t_Baseline * const b )  { baseline = b; }
        //! Returns the pointer to the quantity
        t_Quantity* get_quantity () const { return quantity; }
        //! Returns the pointer to the baseline
        t_Baseline* get_baseline () const { return baseline; }
      
        //! Evaluates the distance between the quantity and the baseline.
        t_Type evaluate();
        //! \brief Evaluates the gradient of the distance between the quantity and the
        //!        baseline.
        virtual void evaluate_gradient(t_Type* const _i_grad);
        //! \brief Returns the gradient in direction \a _pos of the distance
        //!        between the quantity and the baseline 
        virtual t_Type evaluate_one_gradient(types::t_unsigned _pos);
        //! \brief Evaluates the distance between the quantity and the baseline
        //!        and computes its gradient.
        virtual t_Type evaluate_with_gradient(t_Type* const _i_grad);

        //! Initializes the quantity.
        virtual bool init() { return quantity->init(); }
        //! Calls is_taboo() from the quantity.
        virtual bool is_taboo() const  { return quantity->is_taboo(); }

      protected:
        //! \brief Debug: checks that the quantity and the baseline pointers are not
        //!        null.
        void check_pointers(const char *stream) const;
    };

    template<typename T_QUANTITY, typename T_BASELINE>
      Fitness<T_QUANTITY, T_BASELINE> :: Fitness(t_Quantity * const q, t_Baseline * const b) 
      {
        variables = q->get_variables();
        quantity = q; 
        baseline = b; 
        check_pointers("Fitness(t_Quantity *q, t_Baseline *b)");
      }
    
    template<typename T_QUANTITY, typename T_BASELINE>
      void Fitness<T_QUANTITY, T_BASELINE> :: destroy_variables()
      {
        quantity->destroy_variables();
        variables = NULL;
      }

    template<typename T_QUANTITY, typename T_BASELINE>
      void Fitness<T_QUANTITY, T_BASELINE> :: set_quantity ( t_Quantity * const q )
      {
        quantity = q;
        variables = q->get_variables();
      }

    template<typename T_QUANTITY, typename T_BASELINE>
      inline typename Fitness<T_QUANTITY, T_BASELINE> :: t_Type
        Fitness<T_QUANTITY, T_BASELINE> :: evaluate() 
        { 
          check_pointers("evaluate()");
          return (   quantity->evaluate() 
                   - baseline->evaluate( get_concentration() ) );
        }
    template<typename T_QUANTITY, typename T_BASELINE>
      inline void Fitness<T_QUANTITY, T_BASELINE> :: evaluate_gradient(t_Type* const _i_grad) 
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
    template<typename T_QUANTITY, typename T_BASELINE>
      inline typename Fitness<T_QUANTITY, T_BASELINE> :: t_Type
        Fitness<T_QUANTITY, T_BASELINE> :: evaluate_one_gradient(types::t_unsigned _pos)
        {
          t_Type x = get_concentration();
          t_Type value  = quantity->evaluate_one_gradient(_pos);
          t_Type N = t_Type(variables->size());
          t_Type g_shift = baseline->evaluate_gradient( x ) / N;
          return value - g_shift;
        }
    template<typename T_QUANTITY, typename T_BASELINE>
      inline typename Fitness<T_QUANTITY, T_BASELINE> :: t_Type
        Fitness<T_QUANTITY, T_BASELINE> :: evaluate_with_gradient(t_Type* const _i_grad)
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

    template<typename T_QUANTITY, typename T_BASELINE>
      inline void Fitness<T_QUANTITY, T_BASELINE> :: check_pointers(const char *stream) const
      {
        _DODEBUGCODE(
          if ( !quantity or !baseline )
          {
            std::cerr << "Pointers not set in Fitness::" << stream
                      << std::endl;
            exit(0);
          } 
        )
      }

  } // namespace opt
} // namespace LaDa

#endif // _FITNESS_H_
