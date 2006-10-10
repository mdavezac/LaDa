#ifndef _LADA_H_
#define _LADA_H_

//  MotU: Master of the Universe class
//  controls everything, does everything, will cook your breakfast

#include "individual.h"

#include <lamarck/functional_builder.h>
#include <opt/fitness_function.h>
#include <lamarck/ch_template.h>

#include <opt/opt_minimize.h>

#undef min // idiots
#undef max
#include <eo/eoBreed.h>
#include <eo/eoGenOp.h>
#include <eo/eoGenContinue.h>
#include <eo/eoOpContainer.h>
#include <eo/eoReduceMerge.h>
#include <eo/eoDetTournamentSelect.h>
#include <eo/eoGeneralBreeder.h>
#include <eo/utils/eoUpdater.h>
#include <eo/utils/eoState.h>
#include <eo/utils/eoCheckPoint.h>
#undef min // idiots
#undef max

#include "darwin.h"

namespace LaDa
{
  typedef LaDa::Individual<> t_individual ; 
  typedef opt::Fitness_Function<VA_CE::t_functional, VA_CE::CH_Template> t_fitness;
  class MotU : public VA_CE :: Functional_Builder
  {
    protected:
      typedef t_fitness FITNESS;
      const static unsigned LAMARCK;
      const static unsigned DARWIN;
      const static unsigned DEBUG;

    public:
      struct GA // stores all GA parameters
      {
        const static unsigned NO_MINIMIZER;
        const static unsigned WANG_MINIMIZER;
        const static unsigned PHYSICAL_MINIMIZER;
        const static unsigned LINEAR_MINIMIZER;
        const static unsigned SA_MINIMIZER; // simulated annealing at zero T

        double crossover_value; 
        double mutation_value;
        bool sequential_op;
        unsigned tournament_size;
        unsigned pop_size;
        double replacement_rate;
        unsigned max_generations;
        unsigned method;
        bool utter_random;
        bool multistart;
        bool evolve_from_start;
        unsigned minimizer;
        unsigned max_calls;
        double minimize_best;
        bool minimize_offsprings;
        unsigned minimize_best_every;

        // CH stuff
        bool is_one_point_hull;

        // eo stuff
        eoState eostates;
        eoIncrementorParam<unsigned> *nb_generations;
        Darwin<t_individual> *algorithm;

        GA();
        bool Load( TiXmlElement *element );
      };


    private:
      Ising_CE::Structure structure;
      FUNCTIONAL functional;
      FITNESS fitness;
      VA_CE::CH_Template *convex_hull;
      std::string filename;
      std::string xml_filename;
      std::string xmgrace_filename;
      GA ga_params;
      eoPop<t_individual> population;
      unsigned EvalCounter;
      opt::Minimize_Base<FITNESS> *minimizer;


    public:
      MotU() : Functional_Builder(), convex_hull(NULL), filename("input.xml"),
               ga_params(), population(), EvalCounter(0), minimizer(NULL) {}
      bool Load(const std::string &_filename);
      virtual ~MotU(){};
      void run();
      void print_xml();
      void print_xmgrace();

      double evaluate( t_individual &individual ); // polynomial+CS only
      double minimize( t_individual &individual ); // polynomial+CS only
      double evaluate( const double &_x ) // convex hull only
        { return convex_hull->evaluate(_x); }
    
    protected:
      virtual bool Load( TiXmlHandle &handle );
      bool read_CH();
      void run_debug();


      void populate();
      eoCheckPoint<t_individual>* make_checkpoint();
      eoGenOp<t_individual>* make_GenOp();
      eoGenOp<t_individual>* make_recurrent_op(const TiXmlElement * const el);
      eoBreed<t_individual>* make_breeder();
      eoReplacement<t_individual>* make_replacement();
      void make_algo();

      void write_xmgrace_header( std::ofstream &_f);
      void init_convex_hull();
      
  };

} // namespace LaDa

#endif // _LADA_H_
