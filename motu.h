#ifndef _LADA_H_
#define _LADA_H_

//  MotU: Master of the Universe class
//  controls everything, does everything, will cook your breakfast

#include "individual.h"

#include <lamarck/functional_builder.h>
#include <opt/fitness_function.h>
#ifdef ONE_POINT
  #include "one_point_hull.h"
  #define CONVEX_HULL One_Point_Hull
#else 
  #include <lamarck/convex_hull.h>
  #define CONVEX_HULL VA_CE::Convex_Hull
#endif

#define LINEAR_SOLVE
#include <opt/opt_minimize.h>

#include <eo/eo>


namespace LaDa
{
  class MotU : public VA_CE :: Functional_Builder
  {
    public:
      typedef Individual<FUNCTIONAL> t_individual;

    protected:
      typedef opt::Fitness_Function<FUNCTIONAL, CONVEX_HULL> FITNESS;
      const static double ZERO_TOLERANCE;
      const static int LAMARCKIAN_EVOLUTION;
      const static int LAMARCKIAN;
      const static int DARWINISTIC;
      const static int MULTISTART;

      struct GA_Params
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

        GA_Params();
        bool Load( TiXmlElement *element );
      };

    private:
      Ising_CE::Structure structure;
      FUNCTIONAL functional;
      FITNESS fitness;
      opt::Minimize<FITNESS> minimizer;
      CONVEX_HULL convex_hull;
      std::string filename;
      std::string xml_filename;
      std::string xmgrace_filename;
      eoPop<t_individual> population;
      GA_Params ga_params;
      eoState eostates;
      eoIncrementorParam<unsigned> *nb_generations;
      unsigned EvalCounter;
      eoEasyEA<t_individual> *algorithm;


    public:
      MotU() { filename="input.xml"; EvalCounter = 0; };
      MotU(const std::string &_filename);
      ~MotU(){};
      void run();
      void print_xml();
      void print_xmgrace();

      double evaluate( t_individual &individual ); // polynomial only
      double evaluate( const double &_x ) // convex hull only
        { return convex_hull.evaluate(_x); }
    
    protected:
      bool Load( TiXmlHandle &handle );
      bool read_CH();


      void populate();
      eoCheckPoint<t_individual>* make_checkpoint();
      eoGenOp<t_individual>* make_GenOp();
      eoBreed<t_individual>* make_breeder();
      eoReplacement<t_individual>* make_replacement();
      void make_algo();
      
  };

} // namespace LaDa

#endif // _LADA_H_
