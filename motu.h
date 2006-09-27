#ifndef _LADA_H_
#define _LADA_H_

//  MotU: Master of the Universe class
//  controls everything, does everything, will cook your breakfast

#include "individual.h"

#include <lamarck/functional_builder.h>
#include <opt/fitness_function.h>
#include <lamarck/convex_hull.h>

#include <opt/opt_minimize.h>

#include <eo/eo>


namespace LaDa
{
  class MotU : public VA_CE :: Functional_Builder
  {
    public:
      typedef Individual<> t_individual;

    protected:
      typedef opt::Fitness_Function<FUNCTIONAL, VA_CE::Convex_Hull> FITNESS;
      const static unsigned LAMARCK;
      const static unsigned DARWIN;
      const static unsigned DEBUG;
      const static unsigned NO_MINIMIZER;
      const static unsigned WANG_MINIMIZER;
      const static unsigned PHYSICAL_MINIMIZER;
      const static unsigned LINEAR_MINIMIZER;

      struct GA // stores all GA parameters
      {
        double crossover_vs_mutation,
               crossover_probability, 
               mutation_probability;
        bool sequential_op;
        unsigned tournament_size;
        unsigned pop_size;
        double replacement_rate;
        unsigned max_generations;
        unsigned method;
        bool utter_random;
        bool multistart;
        bool evolve_from_start;

        // eo stuff
        eoState eostates;
        eoIncrementorParam<unsigned> *nb_generations;
        eoEasyEA<MotU::t_individual> *algorithm;

        GA();
        bool Load( TiXmlElement *element );
      };


    private:
      Ising_CE::Structure structure;
      FUNCTIONAL functional;
      FITNESS fitness;
      VA_CE::Convex_Hull convex_hull;
      std::string filename;
      std::string xml_filename;
      std::string xmgrace_filename;
      GA ga_params;
      eoPop<t_individual> population;
      unsigned EvalCounter;
      int job_type;
      opt::Minimize_Base<FITNESS> *minimizer;
      unsigned minimizer_type;


    public:
      MotU() : Functional_Builder(), convex_hull(), filename("input.xml"),
               ga_params(), population(), EvalCounter(0)
        { minimizer=NULL; minimizer_type=0; }
      MotU(const std::string &_filename);
      virtual ~MotU(){};
      void run();
      void print_xml();
      void print_xmgrace();

      double evaluate( t_individual &individual ); // polynomial only
      double evaluate( const double &_x ) // convex hull only
        { return convex_hull.evaluate(_x); }
    
    protected:
      bool Load( TiXmlHandle &handle );
      bool read_CH();
      void run_debug();


      void populate();
      eoCheckPoint<t_individual>* make_checkpoint();
      eoGenOp<t_individual>* make_GenOp();
      eoBreed<t_individual>* make_breeder();
      eoReplacement<t_individual>* make_replacement();
      void make_algo();

      void write_xmgrace_header( std::ofstream &_f);
      
  };

} // namespace LaDa

#endif // _LADA_H_
