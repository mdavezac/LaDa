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
  typedef LaDa::Individual<> t_individual ; 
  typedef opt::Fitness_Function<VA_CE::t_functional, VA_CE::CH_Template> t_fitness;
  typedef VA_CE :: CH_Template t_convex_hull;
  class Lamarck : public VA_CE :: Functional_Builder
  {
    public: 
      unsigned EvalCounter;

    protected:
      typedef t_fitness FITNESS;

    private:
      Ising_CE::Structure structure;
      FUNCTIONAL functional;
      t_fitness fitness;
      t_convex_hull *convex_hull;
      std::string filename;
      std::string xml_filename;
      bool is_one_point_hull;
      std::vector< opt::Minimize_Base<t_fitness> *> minimizers;


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
      double evaluate( t_individual & _individual );
      bool minimize( const t_individual & _individual, const unsigned &_minimizer );
    
      t_fitness& get_functional (const t_individual &_individual) 
        { return fitness; }
      
      t_convex_hull& get_convex_hull ()
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

    protected:
      virtual bool Load( TiXmlHandle &handle );
      bool read_CH();
      void run_debug();

      void init_convex_hull();
  };

} // namespace LaDa

#endif // _LADA_H_
