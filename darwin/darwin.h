//
//  Version: $Id$
//
#ifndef  _DARWIN_H_
#define  _DARWIN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <list>

#include <eo/eoPop.h>
#include <eo/eoGenOp.h>
#include <eo/utils/eoState.h>
#include <eo/eoReplacement.h>

#include <tinyxml/tinyxml.h>

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "checkpoints.h"
#include "taboos.h"
#include "objective.h"
#include "store.h"
#include "evaluation.h"
#include "breeder.h"
#include "gatraits.h"

namespace darwin
{
  template< class T_EVALUATOR, class T_GATRAITS = Traits::GA<T_EVALUATOR> >
  class Darwin
  {    
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GATRAITS t_GA_Traits;
    private:
      typedef Darwin<t_Evaluator, t_GA_Traits> t_Darwin;
      typedef typename t_GA_Traits :: t_Individual       t_Individual;
      typedef typename t_GA_Traits :: t_QuantityTraits   t_QuantityTraits;
      typedef typename t_GA_Traits :: t_IndivTraits      t_IndivTraits;
      typedef typename t_IndivTraits :: t_Object         t_Object;
      typedef typename t_IndivTraits :: t_Population     t_Population;
      typedef typename t_IndivTraits :: t_Islands        t_Islands;
      typedef typename t_QuantityTraits :: t_ScalarQuantity   t_ScalarQuantity;
      typedef typename Objective :: Types < t_Evaluator, t_GA_Traits >  t_Objective;
      typedef typename Store :: Types< t_Evaluator, t_GA_Traits >  t_Store;

    protected:
      const static types::t_unsigned SAVE_RESULTS;
      const static types::t_unsigned SAVE_HISTORY;
      const static types::t_unsigned SAVE_POPULATION;

    protected:
      std::string filename;
      std::string evaluator_filename;
      std::string restart_filename;
      std::string save_filename;
      types::t_unsigned tournament_size;
      types::t_unsigned pop_size;
      types::t_unsigned max_generations;
      types::t_unsigned nb_islands;
      types::t_unsigned restarts;
      types::t_unsigned do_save;
      types::t_unsigned do_restart;
      types::t_real     replacement_rate;
      bool do_print_each_call;
      enum { RANDOM_POPULATE, PARTITION_POPULATE } populate_style;

      IslandsContinuator<t_Individual,t_IndivTraits>*          continuator;
      History<t_Individual>*                                   history;
      Taboo_Base<t_Individual>*                                taboos;
      eoGenOp<t_Individual>*                                   breeder_ops;
      Breeder<t_Individual, t_IndivTraits>*                    breeder;
      eoReplacement<t_Individual>*                             replacement;
      typename t_Objective::Vector*                            objective;
      typename t_Store :: Base*                                store;
      Evaluation::Base<t_Evaluator, t_GA_Traits>*              evaluation;

      t_Evaluator   evaluator;
      t_Islands     islands;
      t_Population  offsprings;
      
      eoState eostates;

    public:
      Darwin () : filename("input.xml"), tournament_size(2), pop_size(100),
                  max_generations(0), nb_islands(1), restarts(0), do_save(SAVE_RESULTS),
                  do_restart(0), replacement_rate(0.1), do_print_each_call(false),
                  populate_style(RANDOM_POPULATE), continuator(NULL), history(NULL),
                  taboos(NULL), breeder_ops(NULL), breeder(NULL), replacement(NULL),
                  objective(NULL), store(NULL), evaluation(NULL) {}
      virtual ~Darwin ();

      bool Load( const std::string &_filename );
      void run();

    protected: 
      bool Load_Parameters( const TiXmlElement &_parent );
      void make_History( const TiXmlElement &_parent );
      bool Load_Mating( const TiXmlElement &_parent );
      void Load_Method( const TiXmlElement &_parent );
      void Load_Taboos( const TiXmlElement &_node );
      void Load_CheckPoints (const TiXmlElement &_parent);
      bool Restart(const TiXmlElement &_node);
      bool Restart();
      bool Save(const TiXmlElement &_node);
      bool Save();
      eoGenOp<t_Individual>* make_genetic_op( const TiXmlElement &el,
                                              eoGenOp<t_Individual> *current_op = NULL);
      void make_breeder();
      eoReplacement<t_Individual>* make_replacement();

      void populate ();
#ifdef _MPI
      bool broadcast_islands( mpi::BroadCast &_bc );
      void LoadAllInputFiles(std::string &_input, 
                             std::string &_restart, 
                             std::string &_evaluator );
#endif
      void random_populate ( t_Population &_pop, types::t_unsigned _size);
      void partition_populate ( t_Population &_pop, types::t_unsigned _size);
  };

  template< class T_EVALUATOR, class T_GATRAITS >
    const types::t_unsigned Darwin<T_EVALUATOR,T_GATRAITS> :: SAVE_RESULTS    = 1;
  template< class T_EVALUATOR, class T_GATRAITS >
    const types::t_unsigned Darwin<T_EVALUATOR, T_GATRAITS> :: SAVE_POPULATION = 2;
  template< class T_EVALUATOR, class T_GATRAITS >
    const types::t_unsigned Darwin<T_EVALUATOR, T_GATRAITS> :: SAVE_HISTORY     = 4;

} // namespace darwin

#include "darwin.impl.h"

#endif // _DARWIN_H_

