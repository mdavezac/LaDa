//
//  Version: $Id$
//
#ifndef _MULTIOB_EVALUATION_H_
#define _MULTIOB_EVALUATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list> 

#include <eo/eoPop.h>

#include <tinyxml/tinyxml.h>

#include <mpi/mpi_object.h>
#include <print/stdout.h>

#include <typeinfo>

#include "objective.h"
#include "store.h"
#include "taboos.h"
#include "gatraits.h"

namespace GA
{
/** \ingroup Genetic 
 * @{*/
//! \brief Contains wrapper classes which make call to an Evaluator's functional routines
//! \details The classes in this namespace form the funnel through which all
//! requests from %GA to evaluate any object are sent. The goal is to capture any
//! such request for storage puposes <em>via</em> the classes in the Store
//! namespace, and for taboo (history) purposes <em>via</em> Taboo::History. It
//! also allows us to separate the functional returned by the functional and
//! stored in Base::t_Individual::quantities() from the fitness of an individual,
//! implemented in Fitness and kept in Individual %type objects. Fitness and
//! quantities are linked <em>via</em> the classes implemented in Objective.
//! At this point, there are two classes:
//!  - Base: declares a simple base class with a complete set of storage behaviors
//!  - WithHistory: can also keep tab on all evaluated objects past and present
//! through a Taboo::History object. 
//!  .
//! \sa Store, Store::Base, Taboo::History
namespace Evaluation 
{ 
  //! \brief Funnel through which all calls from %GA to the %Evaluator's
  //! functional routines are routed.
  //! \details Abstract base class.
template< class T_POPULATION >
  class Abstract
  {
    public:
      //! Population type
      typedef T_POPULATION                  t_Population;

    protected:
      //! Individual type
      typedef typename t_Population :: value_type                  t_Individual;
      //! \brief quantity traits pertaining to Virtual Atom minimization
      typedef typename t_Individual :: t_IndivTraits :: t_VA_Traits   t_VATraits;
      //! \brief Gradients type for Virtual Atom minimization
      typedef typename t_VATraits :: t_QuantityGradients           t_QuantityGradients;
      //! \brief Type of the variables used in Virtual Atom minimization 
      //! \details (what is the type of \f$x\f$ in \f$f(x)\f$, where \f$f\f$ is
      //! this evaluation object.
      typedef typename t_VATraits :: t_Type                        t_VA_Type;
      //! \brief Scalar Fitness type
      typedef typename t_Individual :: t_Fitness :: t_ScalarFitness  t_Fitness;
      //! \brief Quantity of the Scalar Fitness 
      typedef typename t_Fitness :: t_Quantity                     t_FitnessQuantity;

    public:
      //! records the number of effective calls to the functional
      types::t_unsigned nb_eval; 
      //! records the number of effective calls to the gradient of the functional
      types::t_unsigned nb_grad; 

      //! Constructor and Initializer
      Abstract () : nb_eval(0), nb_grad(0) {}
      //! Constructor and Initializer
      Abstract ( const Abstract & _c ) : nb_eval( _c.nb_eval ), nb_grad( _c.nb_grad ) {}
      //! Destructor
      virtual ~Abstract() {};

      //! Returns classname.
      virtual std::string className() const { return "Evaluation::Abstract"; }

      //!  calls on the functional to evaluate \a _indiv
      virtual t_FitnessQuantity evaluate( t_Individual &_indiv ) = 0;
      //! calls on the functional to evaluate the (VA) gradients of _indiv
      virtual void evaluate_gradient( t_Individual &_indiv,
                                      t_QuantityGradients& _grad,
                                      t_VA_Type *_i_grad ) = 0;
      //! combines Evaluation::Base::evaluate(t_Individual&) with 
      virtual t_FitnessQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                        t_QuantityGradients& _grad,
                                                        t_VA_Type *_i_grad ) = 0;
      //! calls on the functional to evaluate the (VA) gradients of _indiv
      virtual t_VA_Type evaluate_one_gradient( t_Individual &_indiv,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _pos ) = 0;
      //! Evaluates current populations and offspring populations. 
      virtual void operator()( t_Population &_pop, t_Population &_offspring ) = 0;
      //! Initializations prior to evaluation proper.
      virtual void init( t_Individual &_indiv) = 0;

    protected:
      //! Evaluates a single population.
      virtual void evaluate( t_Population &_pop ) = 0;
  };

  //! \brief Funnel through which all calls from %GA to the %Evaluator's
  //! functional routines are routed.
  //! \details It implements obvious routines such as Base::evaluate(), and
  //! Base::evaluate_gradient, but also a less obvious though necessary storage
  //! ability <em>via</em> through a referenced to a Store::Base objet.
template< class T_GATRAITS >
  class Base : public Abstract< typename T_GATRAITS :: t_Population >
  {
    public:
      //! Should contain all GA definitions \sa Traits::GA
      typedef T_GATRAITS  t_GATraits;

    protected:
      //! Type of the base class.
      typedef Abstract< typename T_GATRAITS :: t_Population > t_Base;
      //! Individual type
      typedef typename t_GATraits :: t_Individual                  t_Individual;
      //! Evaluator type
      typedef typename t_GATraits :: t_Evaluator                   t_Evaluator;
      //! Population type
      typedef typename t_GATraits :: t_Population                  t_Population;
      //! \brief Objective type
      typedef typename Objective :: Types < t_GATraits > :: t_Vector t_Objective;
      //! \brief Store type
      typedef Store :: Base<t_GATraits>                            t_Store;
      //! \brief quantity traits pertaining to Virtual Atom minimization
      typedef typename t_GATraits :: t_VA_Traits                   t_VATraits;
      //! \brief Gradients type for Virtual Atom minimization
      typedef typename t_VATraits :: t_QuantityGradients           t_QuantityGradients;
      //! \brief Type of the variables used in Virtual Atom minimization 
      //! \details (what is the type of \f$x\f$ in \f$f(x)\f$, where \f$f\f$ is
      //! this evaluation object.
      typedef typename t_VATraits :: t_Type                        t_VA_Type;
      //! \brief Scalar Fitness type
      typedef typename t_GATraits :: t_Fitness :: t_ScalarFitness  t_Fitness;
      //! \brief Quantity of the Scalar Fitness 
      typedef typename t_Fitness :: t_Quantity                     t_FitnessQuantity;

    protected:
      //! \brief Evaluator to which all functional calls should be tunneled.
      //! \sa GA::Evaluator, SingleSite::Evaluator, TwoSite::Evaluator,
      //! BandGap::Evaluator, CE::Evaluator, Molecularity::Evaluator
      t_Evaluator *evaluator;
      //! Links quantities and Fitness. \sa Fitness, Objective
      t_Objective *objective; 
      //! Stores results!! \sa Store
      t_Store *store;

    public:
      //! \brief Constructor and Initializer
      Base () : t_Base(), evaluator(NULL), objective(NULL), store(NULL) {}
      //! \brief Constructor and Initializer
      Base   ( t_Evaluator *_eval, t_Objective *_obj, t_Store *_store )
           : t_Base(), evaluator(_eval), objective(_obj), store(_store) {}
      //! \brief Copy Constructor 
      Base   ( const Base<t_Evaluator> &_x )
           : t_Base( _x ), evaluator(_x.evaluator),
             objective(_x.objective), store(_x.store) {}
      //! \brief Destructor
      virtual ~Base() {};

      //! Returns classname.
      virtual std::string className() const { return "Evaluation::Base"; }

    public:
      //! Sets the pointer to the functional interface.
      void set( t_Evaluator &_eval ) { evaluator = &_eval; }
      //! Sets the pointer to the objective interface.
      void set( t_Objective *_obj ) { objective = _obj; }
      //! Sets the pointer to the storage interface.
      void set( t_Store *_stor ) { store = _stor; }
      //! \brief  calls on the functional to evaluate \a _indiv
      //¡ \details All individuals for which the quantities are unknowned
      //! should be evaluated <em>via</em>  this function. The evaluation call
      //! itself is redirected to the Evaluator if and only if the individual is
      //! invalid (meaning its quantities are unknown). The results (e. g. the
      //! quantities) are themselves directed as input to Base::objective. The
      //! result of that is stored into \a _indiv's Fitness object. Finally,
      //! the function also m akes sure that worthy individuals are kept as
      //! results in the object Base::store.
      //! \sa GA::Evaluator::evaluate, SingleSite::Evaluator::evaluate,
      //! TwoSites::Evaluator::evaluate
      virtual t_FitnessQuantity evaluate( t_Individual &_indiv );
      //! \brief calls on the functional to evaluate the (VA) gradients of _indiv
      //! \details Simply links together Objective and the gradients of
      //! Individual::Base::quantities(). Also keeps a tab on the number of
      //! evaluated gradients.
      virtual void evaluate_gradient( t_Individual &_indiv,
                                      t_QuantityGradients& _grad,
                                      t_VA_Type *_i_grad );
      //! \brief combines Evaluation::Base::evaluate(t_Individual&) with 
      //! Evaluation::Base::evaluate_gradient().
      virtual t_FitnessQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                        t_QuantityGradients& _grad,
                                                        t_VA_Type *_i_grad );
      //! \brief calls on the functional to evaluate the (VA) gradients of _indiv
      //! \details Simply links together Objective and the gradients of
      //! Individual::Base::quantities(). Also keeps a tab on the number of
      //! evaluated gradients. Only one gradient is evaluated: \a _pos 
      //! \param _indiv Inidividual from which to get gradient
      //! \param _grad Storage for all gradients, as filled by GA::Evaluator derived classes.
      //! \param _pos index of the gradient to evaluate (eg. position in the
      //! presumably vectorial container Traits::Indiv::t_Object::Container)
      virtual t_VA_Type evaluate_one_gradient( t_Individual &_indiv,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _pos );

    protected:
      //! \brief Wrapper function which calls Base::evaluate(t_Individual &)
      //! for each individual in \a _pop
      //! \details Additionally, this function first calls Base::init to make
      //! sure that Objective::Base::current_individual and
      //! GA::Evaluator::init are called prior to evaluation. In the _MPI
      //! flavor. it also synchronizes Base::store across all processors.
      virtual void evaluate( t_Population &_pop );
    public:
      //! \brief Always evaluates \a _offspring <em>via</em>
      //! Evaluation::Base::evaluate(t_Population&), conditional evaluation of
      //! \a _pop. 
      //! \details \a _pop is assessed only if the objectives are invalid for
      //! some reason. Probably useless since the introduction of Scaling.
      virtual void operator()( t_Population &_pop, t_Population &_offspring );

      //! \brief Correlates Evaluation::Base::objective and
      //! Evaluation::Base::evaluator with \a _indiv
      //! \details This function makes sure that both the Objectives and the
      //! Evaluators know which Individual \a _indiv we are talking about
      void init( t_Individual &_indiv);
  };


  //! \brief Adds History tracking to Evaluation::Base
  //! \details History is managed <em>via</em> the pointer
  //! Evaluation::WithHistory::history to a GA::History object.
  //! History allows us to evaluate any object only once. 
  //! If Evaluation::WithHistory::history is \c NULL then it is ignored.
template< class T_GATRAITS >
  class WithHistory : public Base<T_GATRAITS>
  {
    public:
      //! Should contain all GA definitions \sa Traits::GA
      typedef T_GATRAITS t_GATraits;

    private:
      //! Base class type
      typedef Base<t_GATraits>                                     t_Base;
      //! Individual type
      typedef typename t_GATraits::t_Individual                    t_Individual;
      //! Evaluator type
      typedef typename t_GATraits :: t_Evaluator                   t_Evaluator;
      //! Population type
      typedef typename t_GATraits :: t_Population                  t_Population;
      //! Objective type
      typedef typename Objective :: Types < t_GATraits > :: t_Vector t_Objective;
      //! \brief quantity traits pertaining to Virtual Atom minimization
      typedef typename t_GATraits :: t_VA_Traits                   t_VATraits;
      //! \brief Gradients type for Virtual Atom minimization
      typedef typename t_VATraits :: t_QuantityGradients           t_QuantityGradients;
      //! \brief Type of the variables used in Virtual Atom minimization 
      //! \details (what is the type of \f$x\f$ in \f$f(x)\f$, where \f$f\f$ is
      //! this evaluation object.
      typedef typename t_VATraits :: t_Type                        t_VA_Type;
      //! \brief Store object type
      typedef Store :: Base<t_GATraits>                            t_Store;
      //! \brief History object type
      typedef GA::History<t_Individual>                            t_History;
      //! Scalar Fitness type
      typedef typename t_GATraits :: t_Fitness :: t_ScalarFitness  t_Fitness;
      //! Quantity type of the Scalar Fitness type
      typedef typename t_Fitness :: t_Quantity                     t_FitnessQuantity;

    protected:
      using t_Base :: evaluator;
      using t_Base :: objective;
      using t_Base :: store;
      using t_Base :: nb_eval;
      using t_Base :: nb_grad;
      //! \brief Keeps tab on previously assessed inviduals
      //! \details If \c NULL, then no history tracking is done, in which case
      //! Base and WithHistory are equivalent.
      t_History *history;


    public:
      using t_Base::evaluate;
      using t_Base :: set;

    public:
      //! \brief Constructor.
      WithHistory () : Base<T_GATRAITS>(), history(NULL) {};
      //! \brief Constructor and Initializer
      WithHistory   ( t_Evaluator *_eval, t_Objective *_obj,
                      t_Store *_store, t_History *_hist )
           : Base<T_GATRAITS>(_eval, _obj, _store), history(_hist) {};
      //! \brief Copy Constructor 
      WithHistory   ( const WithHistory<t_Individual> &_c )
                  : Base<T_GATRAITS>( _c ), history(_c.history ) {}
      //! \brief Destructor
      virtual ~WithHistory() {};
    
      //! \brief  calls on the functional to evaluate \a _indiv
      //! \details Performs exactly as Evaluation::Base::evaluate, except that
      //! if Evaluation::WithHistory::history is non-zero, it does keep tracks of
      //! all individuals evaluated in the past. If \a _indiv is the clone of a
      //! previously evaluated individual, its evaluation in
      //! Evaluation::WithHistory::history is used. 
      //! \sa Evaluation::Base::evaluate 
      virtual t_FitnessQuantity evaluate( t_Individual &_indiv );
      //! \brief combines Evaluation::WithHistory::evaluate(t_Individual&) with 
      //! Evaluation::Base::evaluate_gradient().
      virtual t_FitnessQuantity evaluate_with_gradient( t_Individual &_indiv,
                                                        t_QuantityGradients& _grad,
                                                        t_VA_Type *_i_grad );
      //! Sets the pointer to the history interface.
      void set( t_History *_history ) { history = _history; }
  };

} // namespace Evaluation

} // namespace GA

#include "evaluation.impl.h"

/*@}*/

#endif // _RESULTS_H_
