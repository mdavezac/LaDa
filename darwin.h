#ifndef _DARWIN_H
#define _DARWIN_H

//-----------------------------------------------------------------------------

#undef min 
#undef max
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

#include "evaluation.h"
#include "operators.h"
#include "taboo.h"
#include "breeder.h"
#include "checkpoint.h"

namespace LaDa
{

  template<class t_Object, class t_Lamarck>
  class Darwin: public eoAlgo<t_Object>
  {
    protected:
      typedef Darwin<t_Object, t_Lamarck> t_Darwin;

    protected: 
      const static unsigned LAMARCK;
      const static unsigned DARWIN;
      const static unsigned DEBUG;
      const static std::string svn_version;

      // parameters
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
      bool is_one_point_hull;
      std::vector< std::string > print_strings;

      eoState eostates;
      eoIncrementorParam<unsigned> *nb_generations;
      eoPop<t_Object> population, offsprings;

      std::string filename;
      std::string xmgrace_filename;

      eoContinue<t_Object>*          continuator;
      eoPopEvalFunc<t_Object>*       popEval;
      eoBreed<t_Object>*             breeder;
      eoGenOp<t_Object>*             breeder_ops;
      eoReplacement<t_Object>*       replace;
      eoPopAlgo<t_Object>*           extra_popalgo;
      Taboo_Base<t_Object>*          taboos;
      NuclearWinter<t_Object, t_Darwin >* nuclearwinter;

      t_Lamarck *lamarck;

    public:

      Darwin (t_Lamarck *_lam); 
      virtual ~Darwin () {};

      bool Load( const std::string &_filename );

      /// Apply a few generation of evolution to the population.
      virtual void operator()(eoPop<t_Object>& _pop) {}
      void run();

      typename t_Object :: t_Type evaluate( const double &x ) const
        { return lamarck->evaluate( x );  }
      typename t_Object :: t_Type evaluate( t_Object &_object ) const
        { return lamarck->evaluate( _object );  }
      bool minimize( const t_Object &_object, const unsigned &_nb )
        { return lamarck->minimize( _object, _nb);  }
      void print_xmgrace();
      void print_xmgrace( std::string &_str )
        { print_strings.push_back(_str); }
  
      void print_xml()
        { lamarck->print_xml(); };

    protected:
      eoGenOp<t_Object>* make_GenOp(const TiXmlElement *element, std::ofstream &_f);
      eoGenOp<t_Object>* make_genetic_op(const TiXmlElement &el,
                                         std::ofstream &_f,
                                         std::string &_special,
                                         std::string &_base,
                                         eoGenOp<t_Object> *current_op);
      eoBreed<t_Object>* make_breeder();
      eoCheckPoint<t_Object>* make_checkpoint();
      eoReplacement<t_Object>* make_replacement();
      void make_extra_algo();
      void make_algo();
      eoMonOp<t_Object>* make_MonOp(const TiXmlElement &el,
                                    std::ofstream &_f,
                                    std::string &_special,
                                    std::string &_base);
       
      void populate ();
      MinimizationOp< t_Object, Darwin<t_Object, t_Lamarck> >* Load_Minimizer( const TiXmlElement* el,   
                                                                               std::ofstream &_f ); 

      
      bool Load( TiXmlHandle &handle );
      void write_xmgrace_header();
  };

} // namespace LaDa

#include "darwin.impl.h"
#endif

