#ifndef  _DARWIN_H_
#define  _DARWIN_H_

#include <string>
#include <list>

#include <eo/eoPop.h>
#include <eo/eoGenOp.h>
#include <eo/utils/eoState.h>
#include <eo/eoReplacement.h>

#include <tinyxml/tinyxml.h>

#ifdef _MPI
#include<mpi/mpi_object.h>
#endif

#include "checkpoints.h"
#include "taboos.h"
#include "results.h"
#include "breeder.h"
#include "minimizergenop.h"

namespace darwin
{
  template<class T_INDIVIDUAL, class T_EVALUATOR >
  class Darwin
  {    
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR t_Evaluator;
    private:
      typedef Darwin<t_Individual, t_Evaluator> t_Darwin;
      typedef eoPop<t_Individual>  t_Population;
      typedef std::list< t_Population > t_Islands;
      typedef typename t_Individual :: t_Object t_Object;

    protected:
      const static types::t_unsigned SAVE_RESULTS;
      const static types::t_unsigned SAVE_HISTORY;
      const static types::t_unsigned SAVE_POPULATION;

    protected:
      std::string filename;
      std::string evaluator_filename;
      std::string restart_filename;
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

      IslandsContinuator<t_Individual>*                continuator;
      History<t_Individual, std::list<t_Individual> >* history;
      Taboo_Base<t_Individual>*                        taboos;
      Taboo<t_Individual, std::list<t_Individual> >*   agetaboo;
      eoGenOp<t_Individual>*                           breeder_ops;
      NuclearWinter<t_Individual, t_Population >*      nuclearwinter;
//     Colonize<t_Individual>*                          colonize;
//     PopGrowth<t_Individual>*                         popgrowth;
      Breeder<t_Individual>*                           breeder;
      Results<t_Individual, t_Evaluator>*              results;
      eoReplacement<t_Individual>*                     replacement;

      t_Evaluator                                      evaluator;
      t_Islands                                        islands;
      t_Population                                     offsprings;
      
      // lists for cleanup purposes
      typedef std::list< typename minimizer::Base< typename t_Evaluator::t_Functional >* > t_minimizer_list;
      t_minimizer_list minimizers; 

      eoState eostates;

    public:
      Darwin () : filename("input.xml"), tournament_size(2), pop_size(100),
                  max_generations(0), nb_islands(1), restarts(0), do_save(SAVE_RESULTS),
                  do_restart(0), replacement_rate(0.1), do_print_each_call(false),
                  populate_style(RANDOM_POPULATE), continuator(NULL), history(NULL),
                  taboos(NULL), breeder_ops(NULL), nuclearwinter(NULL), breeder(NULL),
                  results(NULL){};
      virtual ~Darwin ();

      bool Load( const std::string &_filename );
      void run();

    protected: 
      bool Load_Parameters( const TiXmlElement &_parent );
      void Load_History( const TiXmlElement &_parent );
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
      eoReplacement<T_INDIVIDUAL>* make_replacement();

      void populate ();
#ifdef _MPI
      bool broadcast( mpi::BroadCast &_bc );
      bool broadcast_islands( mpi::BroadCast &_bc );
#endif
      void random_populate ( t_Population &_pop, types::t_unsigned _size);
      void partition_populate ( t_Population &_pop, types::t_unsigned _size);
  };

  template< class T_INDIVIDUAL, class T_EVALUATOR >
    const types::t_unsigned Darwin<T_INDIVIDUAL, T_EVALUATOR> :: SAVE_RESULTS    = 1;
  template< class T_INDIVIDUAL, class T_EVALUATOR >
    const types::t_unsigned Darwin<T_INDIVIDUAL, T_EVALUATOR> :: SAVE_POPULATION = 2;
  template< class T_INDIVIDUAL, class T_EVALUATOR >
    const types::t_unsigned Darwin<T_INDIVIDUAL, T_EVALUATOR> :: SAVE_HISTORY     = 4;

} // namespace darwin

#include "darwin.impl.h"

#endif // _DARWIN_H_

