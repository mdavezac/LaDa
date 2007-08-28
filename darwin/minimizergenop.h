//
//  Version: $Id$
//
#ifndef _DARWIN_MINIMIZER_GENOP_H_
#define _DARWIN_MINIMIZER_GENOP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/eoGenOp.h>

#include "opt/opt_function_base.h"
#include "evaluation.h"
#include "objective.h"
#include "gatraits.h"
#include "minimizergenop.h"
#include "opt/va_minimizer.h"

namespace darwin 
{

  template< class T_EVALUATOR, class T_GA_TRAITS > class SaveStateIndividual;
  template< class T_EVALUATOR, class T_GA_TRAITS > class MinimizerGenOp;

  template<class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
  class Minimizer_Functional : public eoFunctorBase, 
                               public function::Base < typename T_GA_TRAITS :: t_VA_Traits :: t_Type,
                                                       typename T_GA_TRAITS :: t_VA_Traits :: t_Container >
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef eoPop<t_Individual> t_Population;
      typedef typename t_GA_Traits :: t_VA_Traits t_VA_Traits;
      typedef typename t_VA_Traits :: t_Type t_VA_Type;
      typedef typename t_VA_Traits :: t_Container t_VA_Container;
      typedef typename function::Base< t_VA_Type, t_VA_Container > t_Base;
      typedef Evaluation::Base<t_Evaluator, t_GA_Traits> t_Evaluation;
      typedef darwin::Taboo_Base<t_Individual> t_Taboo;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_GA_Traits :: t_QuantityTraits t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      typedef SaveStateIndividual<t_Evaluator, t_GA_Traits > t_SaveState;


    protected:
      using t_Base :: variables;

    protected:
      t_Evaluation *evaluation;
      t_Taboo *taboo;
      t_Individual *current_indiv;
      t_QuantityGradients gradients;
    public:
      t_SaveState savestate;
     
    public:
      Minimizer_Functional( t_Evaluation &_r, t_Taboo &_t ) : evaluation(&_r), taboo(&_t) {};
      t_VA_Type evaluate()
      {
        current_indiv->invalidate();
        t_VA_Type result = evaluation->evaluate( *current_indiv ); 
        return result;
      }
      t_VA_Type evaluate_with_gradient( t_VA_Type *_i_grad )
      {
        current_indiv->invalidate();
        return (t_VA_Type) evaluation->evaluate_with_gradient( *current_indiv, gradients, _i_grad ); 
      }
      void evaluate_gradient( t_VA_Type *_i_grad )
        { evaluation->evaluate_gradient( *current_indiv, gradients, _i_grad ); }
      t_VA_Type evaluate_one_gradient( types::t_unsigned _pos )
        { return evaluation->evaluate_one_gradient( *current_indiv, gradients, _pos ); }
      bool is_taboo() const
        { return taboo ? (*taboo)( *current_indiv ): false; }
      bool init( t_Individual & _indiv)
      {
        savestate.init(_indiv);
        evaluation->init(_indiv); 
        variables = &_indiv.Object().Container();
        gradients.resize( variables->size(), _indiv.quantities() );
        Traits::zero_out( gradients );
        current_indiv = &_indiv;
        return true;
      }
      bool init() { return true; }
        
  };

  template< class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR > >
  class SaveStateIndividual
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GA_Traits :: t_QuantityTraits t_QuantityTraits;
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      typedef Minimizer_Functional<t_Evaluator, t_GA_Traits> t_MFunctional;

    protected:
      t_Individual *current_indiv;
      t_Quantity quantity;
      darwin::Fitness fitness;

    public:
      SaveStateIndividual() : current_indiv(NULL) {};
      SaveStateIndividual( t_Individual &_c) : current_indiv(_c) {};
      SaveStateIndividual   ( const SaveStateIndividual &_c) 
                          : current_indiv(_c.current_indiv), quantity( _c.quantity ),
                            fitness( _c.fitness ) {};
      ~SaveStateIndividual() {};
      
      void save()
      {
        if ( not current_indiv ) return;
        quantity = current_indiv->quantities();
        fitness =  current_indiv->fitness();
      }
      void reset() const 
      {
        if ( not current_indiv ) return;
        current_indiv->quantities() = quantity;
        current_indiv->fitness() = fitness;
      }
      void init( t_Individual &_indiv ) { current_indiv = &_indiv; }
  };

  template< class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
  class MinimizerGenOp : public eoGenOp<typename T_GA_TRAITS :: t_Individual>
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;
    protected:
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
      typedef typename t_IndivTraits :: t_VA_Traits t_VA_Traits;
      typedef typename t_VA_Traits :: t_Functional t_Functional;
      typedef typename ::minimizer::Base< t_Functional > t_Minimizer;
      typedef Minimizer_Functional<t_Evaluator, t_GA_Traits> t_MFunctional;
      typedef  SaveStateIndividual<t_Evaluator, t_GA_Traits>  t_SaveState;

    protected:
      t_Minimizer *minimizer;
      t_MFunctional functional;

    public:
      explicit
        MinimizerGenOp   ( t_MFunctional &_r )
                       : minimizer(NULL), functional( _r ) {};
      ~MinimizerGenOp () { if ( minimizer ) delete minimizer; }

      unsigned max_production(void) { return 1; } 
   
      void apply(eoPopulator<t_Individual>& _pop)
      {
        functional.init( *_pop );
        (*minimizer)( functional );
        functional.evaluate();
        (*_pop).invalidate();
      }
      virtual std::string className() const {return "Darwin::MinimizerGenOp";}

      bool Load( const TiXmlElement &_node );
  };
  template< class T_EVALUATOR, class T_GA_TRAITS >
  bool MinimizerGenOp<T_EVALUATOR, T_GA_TRAITS> :: Load( const TiXmlElement &_node )
  {
    std::string name = _node.Value();
    if ( name.compare("Minimizer") != 0 )
      return false;
   
    if ( not _node.Attribute("type") )
      return false;
    
    name = _node.Attribute("type");
   
    if ( name.compare("VA") == 0 )
    {
      darwin::printxmg << darwin::PrintXmg::comment << "VA optimizer" << darwin::PrintXmg::endl;
      // pointer is owned by caller !!
      // don't deallocate
      minimizer =  new ::minimizer::VA<t_Functional, t_SaveState>( _node, functional.savestate );
    }
    else if ( name.compare("SA") == 0 )
    {
      darwin::printxmg << darwin::PrintXmg::comment << "SA optimizer" << darwin::PrintXmg::endl;
      // pointer is owned by caller !!
      // don't deallocate
      minimizer = new ::minimizer::VA<t_Functional, t_SaveState>( _node, functional.savestate );
    }
    else if ( name.compare("Beratan") == 0 )
    {
      darwin::printxmg << darwin::PrintXmg::comment << "Beratan optimizer" << darwin::PrintXmg::endl;
      // pointer is owned by caller !!
      // don't deallocate
      minimizer = new ::minimizer::Beratan<t_Functional>( _node );
    }
   
    return minimizer != NULL;
   
  }


} // namespace darwin

#endif // _DARWIN_MINMIZER_GENOP_H_
