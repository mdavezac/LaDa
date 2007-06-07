#ifndef _DARWIN_H
#define _DARWIN_H

//-----------------------------------------------------------------------------

#include <complex>

#include <eo/apply.h>
#include <eo/eoAlgo.h>
#include <eo/eoPopAlgo.h>
#include <eo/eoPopEvalFunc.h>
#include <eo/eoContinue.h>
#include <eo/eoSelect.h>
#include <eo/eoTransform.h>
#include <eo/eoBreed.h>
#include <eo/eoMergeReduce.h>
#include <eo/eoReplacement.h>
#include <eo/eoGenOp.h>
#include <eo/eoGenContinue.h>
#include <eo/utils/eoUpdater.h>
#include <eo/utils/eoState.h>
#include <eo/utils/eoCheckPoint.h>

#include <tinyxml/tinyxml.h>

#include <list>

#include "evaluation.h"
#include "operators.h"
#include "taboo.h"
#include "breeder.h"
#include "checkpoint.h"
#include "colonize.h"
#include "popgrowth.h"

#include "eotypes.h"

#include <opt/types.h>
using namespace types;

namespace LaDa
{

  template<class T_OBJECT, class T_LAMARCK>
  class Darwin: public eoAlgo<T_OBJECT>
  {
    public:
      typedef T_OBJECT t_Object;
      typedef T_LAMARCK t_Lamarck;
      typedef Darwin<t_Object, t_Lamarck> t_Darwin;
      typedef typename std::list< eoPop<t_Object> > t_Islands;
      typedef typename t_Lamarck :: t_GA_Functional t_Functional;

    protected: 
      const static t_unsigned LAMARCK;
      const static t_unsigned DARWIN;
      const static t_unsigned DEBUG;
      const static t_unsigned RANDOM_POPULATE;
      const static t_unsigned PARTITION_POPULATE;
      const static std::string svn_version;

      // parameters
      t_real crossover_value; 
      t_real mutation_value;
      t_real replacement_rate;
      t_real minimize_best;
      bool sequential_op;
      bool utter_random;
      bool multistart;
      bool evolve_from_start;
      bool minimize_offsprings;
      bool is_one_point_hull;
      bool do_print_nb_calls;
      eotypes::t_unsigned tournament_size;
      eotypes::t_unsigned pop_size;
      t_unsigned max_generations;
      t_unsigned method;
      t_unsigned minimizer;
      t_unsigned max_calls;
      t_unsigned minimize_best_every;
      t_unsigned populate_style;
      std::vector< std::string > print_strings;

      eoState eostates;
      t_unsigned nb_islands;
      t_Islands islands;
      eoPop<t_Object> offsprings;

      std::string filename;
      std::string xmgrace_filename;

      IslandsContinuator<t_Object>*  continuator;
      eoPopEvalFunc<t_Object>*       popEval;
      Breeder<t_Object>*             breeder;
      eoGenOp<t_Object>*             breeder_ops;
      eoReplacement<t_Object>*       replace;
      eoPopAlgo<t_Object>*           extra_popalgo;
      Taboo_Base<t_Object>*          taboos;
      Taboo<t_Object, std::list<t_Object> > *agetaboo;
      OffspringTaboo<t_Object, std::list<t_Object> > *pathtaboo;
      History<t_Object, std::list<t_Object> > *history;
      NuclearWinter<t_Darwin >* nuclearwinter;
      Colonize<t_Object> *colonize;
      PopGrowth<t_Object> *popgrowth;
      Evaluation<t_Darwin> *evaluation;

      t_Lamarck *lamarck;

    public:

      Darwin (t_Lamarck *_lam); 
      virtual ~Darwin () {};

      bool Load( const std::string &_filename );

      /// Apply a few generation of evolution to the population.
      virtual void operator()(eoPop<t_Object>& _pop) {}
      void run();

      typename t_Object :: t_Type evaluate( const t_real &x ) const
        { return lamarck->evaluate( x );  }
      typename t_Object :: t_Type evaluate( const t_Object &_object ) const
        { return lamarck->evaluate( _object );  }
      bool minimize( const t_Object &_object, const t_unsigned &_nb )
        { return lamarck->minimize( _object, _nb);  }
      void print_xmgrace(bool is_last_call);
      void print_xmgrace( std::string &_str )
        { print_strings.push_back(_str); }
  
      void print_xml()
        { lamarck->print_xml(); };
      
      void fourrier_transform( const t_Object &_object, 
                               std::vector< std::complex< typename t_Object :: t_Type> > &_fourrier )
        { return lamarck->fourrier_transform( _object, _fourrier ); }
      t_Functional &get_functional( const t_Object &_object ) 
        { return lamarck->get_GA_functional( _object ); }
      void add_to_convex_hull( const t_Object &_object ) const
        { lamarck->add_to_convex_hull(_object); }

    protected:
      eoGenOp<t_Object>* make_GenOp(const TiXmlElement *element, std::ofstream &_f);
      eoGenOp<t_Object>* make_genetic_op(const TiXmlElement &el,
                                         std::ofstream &_f,
                                         std::string &_special,
                                         std::string &_base,
                                         eoGenOp<t_Object> *current_op);
      void make_breeder();
      void make_checkpoint();
      void make_taboos();
      void make_breeder_ops();
      eoReplacement<t_Object>* make_replacement();
      void make_extra_algo();
      void make_algo();
      void make_colonize();
      void make_popgrowth();
      void make_history();
       
      void random_populate ();
      void partition_populate ();
      MinimizationOp<t_Object, t_Darwin>* Load_Minimizer( const TiXmlElement* el,   
                                                          std::ofstream &_f ); 

      
      bool Load( TiXmlHandle &handle );
      void write_xmgrace_header();

  };

} // namespace LaDa

#include "darwin.impl.h"
#endif

