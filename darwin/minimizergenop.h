//
//  Version: $Id$
//
#ifndef _DARWIN_MINIMIZER_GENOP_H_
#define _DARWIN_MINIMIZER_GENOP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/eoGenOp.h>

#include <opt/function_base.h>
#include <opt/va_minimizer.h>
#include <print/xmg.h>

#include "evaluation.h"
#include "objective.h"
#include "fitness.h"
#include "gatraits.h"
#include "minimizergenop.h"

namespace GA 
{

  template< class T_GATRAITS > class SaveStateIndividual;
  template< class T_GATRAITS > class MinimizerGenOp;
 

  /** \brief Wrapper to create a functional for minimization around an Evaluation object.
      \details This class redefines the evaluation subroutines with something
               Minimizer objects can deal with. It can also "save" the state of
               an individual for use with taboos in minimizers that accept
               them. Note that it is the scalar fitness which is minimized
               according to this wrapper. */
  template<class T_GATRAITS>
  class Minimizer_Functional : public eoFunctorBase, 
                               public function::Base < typename T_GATRAITS :: t_VA_Traits :: t_Type,
                                                       typename T_GATRAITS :: t_VA_Traits :: t_Container >
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
      //! The base class of evaluation objects.
      typedef Evaluation::Base<t_GATraits>                         t_Evaluation;
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
      //! A pointer to a taboo instance
      t_Taboo *taboo;
      //! A pointer to the current individual being minimized.
      t_Individual *current_indiv;
      //! Gradients for minimization
      t_QuantityGradients gradients;
    public:
      //! Functor for saving the current state.
      t_SaveState savestate;
     
    public:
      //! Constructor and initializer
      Minimizer_Functional( t_Evaluation &_r, t_Taboo &_t ) : evaluation(&_r), taboo(&_t) {};
      //! Evaluates and returns the evaluation
      t_VA_Type evaluate() { return (t_VA_Type) evaluation->evaluate( *current_indiv ); }
      //! \brief Evaluates and returns the evaluation. Also computes the gradient and
      //!        stores in in \a _i_grad.
      //! \pre \a _i_grad should point to a valid and sufficiently large memory block.
      t_VA_Type evaluate_with_gradient( t_VA_Type *_i_grad );
      //! \brief Computes the gradient and stores it in \a _i_grad.
      //! \pre \a _i_grad should point to a valid and sufficiently large memory block.
      void evaluate_gradient( t_VA_Type *_i_grad )
        { evaluation->evaluate_gradient( *current_indiv, gradients, _i_grad ); }
      //! Returns the gradient evaluated in direction \a _pos
      t_VA_Type evaluate_one_gradient( types::t_unsigned _pos )
        { return evaluation->evaluate_one_gradient( *current_indiv, gradients, _pos ); }
      //! Returns true if the current state is taboo
      bool is_taboo() const { return taboo ? (*taboo)( *current_indiv ): false; }
      //! Does all initialization  neeed before minimizing.
      bool init( t_Individual & _indiv);
      //! Invalidates the current invididual so that it can evaluated anew.
      void invalidate() { current_indiv->invalidate(); }
      //! Returns true.
      bool init() { return true; }
        
  };

  //! \brief Functor for saving current minimization state
  //! \details For any number of reasons the minimizer may want to go back
  //!          to a previous state. This functor is used by minimization
  //!          functors to save the current state and possibly reload it later
  //!          on. In practice, and in the case Virtual-Atom like minimization,
  //!          saving the "state" means saving both the quantity and the
  //!          fitness. It is assumed that the variables of the functionals for
  //!          the saved state is recovered by the minizer itself. For an
  //!          example, see minimizer::VA.
  template< class T_GATRAITS >
  class SaveStateIndividual
  {
    public:
      //! All %types relevant to %GA.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! The type of the individual.
      typedef typename t_GATraits :: t_Individual                  t_Individual;
      //! All %types pertaining to the quantity
      typedef typename t_GATraits :: t_QuantityTraits              t_QuantityTraits;
      //! The type of the quantity
      typedef typename t_QuantityTraits :: t_Quantity              t_Quantity;
      //! The type of the scalar quantity
      typedef typename t_QuantityTraits :: t_ScalarQuantity        t_ScalarQuantity;
      //! The type of the scalar fitness.
      typedef typename t_GATraits :: t_Fitness :: t_ScalarFitness  t_Fitness;

    protected:
      //! Pointer to the individual begin  minimized
      t_Individual *current_indiv;
      //! Instance of t_Quantity for saving the state
      t_Quantity quantity;
      //! Instance of t_Fitness for saving the state
      t_Fitness fitness;

    public:
      //! Constructor
      SaveStateIndividual() : current_indiv(NULL) {};
      //! Constructor and Initializer
      SaveStateIndividual( t_Individual &_c) : current_indiv(_c) {};
      //! Copy Constructor
      SaveStateIndividual   ( const SaveStateIndividual &_c) 
                          : current_indiv(_c.current_indiv), quantity( _c.quantity ),
                            fitness( _c.fitness ) {};
      //! Destructor
      ~SaveStateIndividual() {};
      
      //! Saves the state.
      void save();
      //! resets the state to the saved state.
      void reset() const;
      //! sets SaveStateIndividual::current_indiv to \a _indiv.
      void init( t_Individual &_indiv ) { current_indiv = &_indiv; }
  };

  //! \brief Wrapper around a minimizer to create a genetic operator.
  //! \details This operator loads itself from XML. Two types of minimizers are
  //!          accepted: minimizer::VA and minimizer :: Beratan.
  template< class T_GATRAITS >
  class MinimizerGenOp : public eoGenOp<typename T_GATRAITS :: t_Individual>
  {
    public:
      //! All %types relevant to %GA.
      typedef T_GATRAITS                                           t_GATraits;
    protected:
      //! The wrapper functional to minimize
      typedef Minimizer_Functional<t_GATraits>                     t_Functional;
      //! Type of the individual
      typedef typename t_GATraits :: t_Individual                  t_Individual;
      //! All %types relevant to the individual
      typedef typename t_GATraits :: t_IndivTraits :: t_VA_Traits  t_VA_Traits;
      //! All %types relevant to the minimization
//     typedef typename t_VA_Traits :: t_Functional                 t_Functional;
//     //! The base minimizer type
      typedef typename Minimizer::Base< t_Functional >             t_Minimizer;
      //! A functor capable of saving a current state in a minimization functor.
      typedef SaveStateIndividual<t_GATraits >                     t_SaveState;

    protected:
      //! A pointer to a minizer
      t_Minimizer *minimizer;
      //! The wrapper functional
      t_Functional functional;

    public:
      //! Constructor and Initializer
      explicit
        MinimizerGenOp   ( t_Functional &_r )
                       : minimizer(NULL), functional( _r ) {};
      //! Destructor
      ~MinimizerGenOp () { if ( minimizer ) delete minimizer; }

      //! The number of created offspring.
      unsigned max_production(void) { return 1; } 
   
      //! This is a genetic functor
      void apply(eoPopulator<t_Individual>& _pop);

      //! Returns "GA::MinimizerGenOp"
      virtual std::string className() const {return "GA::MinimizerGenOp";}

      //! Loads the minimizer from XML.
      bool Load( const TiXmlElement &_node );
  };

} // namespace GA

#include "minimizergenop.impl.h"

#endif // _DARWIN_MINMIZER_GENOP_H_
