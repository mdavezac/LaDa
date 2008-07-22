#ifndef _LAMARCK_H_
#define _LAMARCK_H_

#include <string>
#include <stdio.h>

#include <lamarck/functional_builder.h>
#include <opt/fitness_function.h>
#ifdef ONE_POINT
  #include "one_point_hull.h"
  #define CONVEX_HULL One_Point_Hull
#else 
  #include <lamarck/convex_hull.h>
  #define CONVEX_HULL VA_CE::Convex_Hull
#endif

// #define WANG_CONSTRAINTS 
#define LINEAR_SOLVE
#include <opt/opt_minimize.h>


#undef min
#undef max
#undef MIN
#undef MAX
#include "atat/xtalutil.h"
#undef min
#undef max
#undef MIN
#undef MAX

#include "solver.h"


class Lamarck : public VA_CE::Functional_Builder
{
  protected: 
    typedef opt::Fitness_Function<FUNCTIONAL, CONVEX_HULL> FITNESS;
    const static double ZERO_TOLERANCE;
    const static int LAMARCKIAN_EVOLUTION;
    const static int LAMARCKIAN;
    const static int DARWINISTIC;
    const static int MULTISTART;

  public:
    const static int XMGRACE_FORMAT;
    const static int XML_FORMAT;


  public:
    bool CalculatingOldPop;
    Ising_CE::Structure structure;

  protected:
    FUNCTIONAL functional;
    FITNESS fitness;
    opt::Minimize<FITNESS> minimizer;
    CONVEX_HULL convex_hull;
    std::string filename;
    std::string xml_filename;
    std::string xmgrace_filename;
    bool centered;
    #ifdef WANG_CONSTRAINTS
      double weight;
    #endif
    int GA_style;

    // parameters
    double random_range;

  public:
    bool convex_hull_has_changed;
    int nb_iterations;

  public:
    Lamarck(){};
    Lamarck(const std::string &_filename);
    ~Lamarck();

    bool Load(const std::string &_filename);
    double evaluate(AtomicConfig * pAtoms);
    double optimize(AtomicConfig * pAtoms);
    double no_optimize(AtomicConfig * pAtoms);
    void iaga_setup(SolverContext* solver);
      
    double evaluate_convex_hull(const double &concentration )
      { return convex_hull.evaluate(concentration); }
    void increment_iteration_count()
      { nb_iterations++; }
    void print_out(int format);
  

  protected:
    bool Load( TiXmlHandle &handle );
    bool read_CH();
    void print_xml();
    void print_xmgrace();

#ifdef WANG_CONSTRAINTS
  protected:
    double find_best_vertex( const AtomicConfig *pAtoms );
#endif

};

#endif
