//
//  Version: $Id: minimizergenop.h 844 2008-11-08 01:22:54Z davezac $
//
#ifndef _LADA_GA_OPERATOR_DETAILS_FUNCTIONALWRAPPER_H_
#define _LADA_GA_OPERATOR_DETAILS_FUNCTIONALWRAPPER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <opt/function_base.h>
#include "savestates.h"

namespace LaDa
{
  namespace GA 
  {

    namespace Operator
    {
      namespace details
      {
        /** \brief Wrapper to create a functional for minimization around an Evaluation object.
            \details This class redefines the evaluation subroutines with something
                     Minimizer objects can deal with. It can also "save" the state of
                     an individual for use with taboos in minimizers that accept
                     them. Note that it is the scalar fitness which is minimized
                     according to this wrapper. */
        template<class T_GATRAITS>
        class FuncWrapper : public function::Base
                            < 
                              typename T_GATRAITS :: t_VA_Traits :: t_Type,
                              typename T_GATRAITS :: t_VA_Traits :: t_Container 
                            >
        {
          public:
            //! All %types relevant to %GA.
            typedef T_GATRAITS t_GATraits;
          protected:
            //! Type of the individual
            typedef typename t_GATraits :: t_Individual                  t_Individual;
            //! %Types pertinent to minimizations.
            typedef typename t_GATraits :: t_VA_Traits                   t_VA_Traits;
            //! The type of scalar component of the minimization argument
            typedef typename t_VA_Traits :: t_Type                       t_VA_Type;
            //! The type the minimization argument
            typedef typename t_VA_Traits :: t_Container                  t_VA_Container;
            //! The base of all functions...
            typedef typename function::Base< t_VA_Type, t_VA_Container > t_Base;
            //! The population container.
            typedef typename t_GATraits :: t_Population                  t_Population;
            //! The base class of evaluation objects.
            typedef Evaluation::Abstract<t_Population>                   t_Evaluation;
            //! The base calss of the taboo objects
            typedef GA::Taboo_Base<t_Individual>                         t_Taboo;
            //! The type of the gradients used in minimization
            typedef typename t_VA_Traits :: t_QuantityGradients          t_QuantityGradients;
            //! All %types relevant to quantity being searched
            typedef typename t_GATraits :: t_QuantityTraits              t_QuantityTraits;
            //! The type of the quantity being searched
            typedef typename t_QuantityTraits :: t_Quantity              t_Quantity;
            //! The type of the \e scalar quantity being searched
            typedef typename t_QuantityTraits :: t_ScalarQuantity        t_ScalarQuantity;
            //! A functor capable of saving a current state in a minimization functor.
            typedef SaveStateIndividual<t_GATraits >                     t_SaveState;
 
 
          protected:
            using t_Base :: variables;
 
          protected:
            //! A pointer to an evaluation instance
            t_Evaluation *evaluation;
            //! A callback to a taboo.
            boost::function<bool(const t_Individual& )> taboo;
            //! A pointer to the current individual being minimized.
            t_Individual *current_indiv;
            //! Gradients for minimization
            t_QuantityGradients gradients;
          public:
            //! Functor for saving the current state.
            t_SaveState savestate;
           
          public:
            //! Constructor and initializer
            Minimizer_Functional   ( t_Evaluation &_r, t_Taboo *_t )
                                 : evaluation(&_r) { set_taboo( _t ); }
            //! Evaluates and returns the evaluation
            t_VA_Type evaluate();
            //! \brief Evaluates and returns the evaluation. Also computes the gradient and
            //!        stores in in \a _i_grad.
            //! \pre \a _i_grad should point to a valid and sufficiently large memory block.
            t_VA_Type evaluate_with_gradient( t_VA_Type *_i_grad );
            //! \brief Computes the gradient and stores it in \a _i_grad.
            //! \pre \a _i_grad should point to a valid and sufficiently large memory block.
            void evaluate_gradient( t_VA_Type *_i_grad )
              { evaluation->evaluate_gradient( *current_indiv, gradients, _i_grad ); }
            //! Returns the gradient evaluated in direction \a _pos
            t_VA_Type evaluate_one_gradient( types::t_unsigned _pos );
            //! Returns true if the current state is taboo
            bool is_taboo() const { return taboo ? (*taboo)( *current_indiv ): false; }
            //! Does all initialization  neeed before minimizing.
            bool init( t_Individual & _indiv);
            //! Invalidates the current invididual so that it can evaluated anew.
            void invalidate() { current_indiv->invalidate(); }
            //! Returns true.
            bool init() { return true; }
              
        };
 
      } // namespace details
    } // namespace Operator

  } // namespace GA
} // namespace LaDa

#include "minimizergenop.impl.h"

#endif // _DARWIN_MINMIZER_GENOP_H_
