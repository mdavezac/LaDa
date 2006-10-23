#ifndef _LADA_H_
#define _LADA_H_

//  MotU: Master of the Universe class
//  controls everything, does everything, will cook your breakfast

#include "individual.h"

#include <lamarck/functional_builder.h>
#include <lamarck/ch_template.h>

#include <opt/fitness_function.h>
#include <opt/opt_minimize.h>


namespace LaDa
{
  class Lamarck : public VA_CE :: Functional_Builder
  {
    public:
      using VA_CE::Functional_Builder::t_VA_Functional;
      typedef LaDa::Individual<> t_Individual; 
      typedef VA_CE :: CH_Template t_Convex_Hull;
      typedef opt::Fitness_Function<t_VA_Functional, t_Convex_Hull> t_GA_Functional;

    public: 
      unsigned EvalCounter;

    private:
      Ising_CE::Structure structure;
      t_VA_Functional functional;
      t_GA_Functional fitness;
      t_Convex_Hull *convex_hull;
      std::string filename;
      std::string xml_filename;
      bool is_one_point_hull;
      std::vector< opt::Minimize_Base<t_GA_Functional> *> minimizers;


    public:
      Lamarck() : Functional_Builder(), EvalCounter(0), convex_hull(NULL), filename("input.xml"),
                  is_one_point_hull(false), minimizers()
      {
        minimizers.reserve(5);
      }
      virtual ~Lamarck();
      bool Load(const std::string &_filename);
      void print_xml();
      void print_xmgrace();

      double evaluate( const double &_x ) // convex hull only
        { return convex_hull->evaluate(_x); }
      double evaluate( t_Individual & _individual );
      bool minimize( const t_Individual & _individual, const unsigned &_minimizer );
    
      t_GA_Functional& get_functional (const t_Individual &_individual) 
        { return fitness; }
      
      t_Convex_Hull& get_convex_hull ()
        { return *convex_hull; }

      unsigned add_minimizer( unsigned type, unsigned n);

      unsigned get_pb_size () const
        { return structure.atoms.size(); }

      void write_xmgrace_header( std::ofstream &_f )
      {
        if ( is_one_point_hull )
          _f << "# One Point Convex Hull " << std::endl;
        else
          _f << "# N-Point Convex Hull " << std::endl;
      }
      void print_xmgrace( std::ofstream &_f, bool _print_ch );

      t_GA_Functional& get_GA_functional( const t_Individual &_object )
        { return fitness; }
      void add_to_convex_hull( const t_Individual &_indiv );

      bool kCrossover( t_Individual  &_offspring, const t_Individual &_parent);
    protected:
      virtual bool Load( TiXmlHandle &handle );
      bool read_CH();
      void run_debug();

      void init_convex_hull();
  };

} // namespace LaDa

#endif // _LADA_H_
