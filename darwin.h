#ifndef _DARWIN_H
#define _DARWIN_H

//-----------------------------------------------------------------------------

#include <apply.h>
#include <eoAlgo.h>
#include <eoPopAlgo.h>
#include <eoPopEvalFunc.h>
#include <eoContinue.h>
#include <eoSelect.h>
#include <eoTransform.h>
#include <eoBreed.h>
#include <eoMergeReduce.h>
#include <eoReplacement.h>
#include <eo/eoGenOp.h>
#include <eo/eoGenContinue.h>
#include <eo/utils/eoUpdater.h>
#include <eo/utils/eoState.h>
#include <eo/utils/eoCheckPoint.h>

#include <tinyxml/tinyxml.h>

#include "evaluation.h"
#include "operators.h"

namespace LaDa
{

  template<class t_Object, class t_Lamarck>
  class Darwin: public eoAlgo<t_Object>
  {
    protected: 
      const static unsigned LAMARCK;
      const static unsigned DARWIN;
      const static unsigned DEBUG;

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

      eoState eostates;
      eoIncrementorParam<unsigned> *nb_generations;
      eoPop<t_Object> population;

      std::string filename;
      std::string xmgrace_filename;

      eoContinue<t_Object>*          continuator;
      eoPopEvalFunc<t_Object>*       popEval;
      eoBreed<t_Object>*             breed;
      eoReplacement<t_Object>*       replace;
      eoPopAlgo<t_Object>*           extra_popalgo;

      t_Lamarck *lamarck;

    public:

      Darwin (t_Lamarck *_lam); 

      bool Load( const std::string &_filename );

      /// Apply a few generation of evolution to the population.
      virtual void operator()(eoPop<t_Object>& _pop) {}
      void run();

      typename t_Object :: TYPE evaluate( const double &x ) const
        { return lamarck->evaluate( x );  }
      typename t_Object :: TYPE evaluate( t_Object &_object ) const
        { return lamarck->evaluate( _object );  }
      bool minimize( const t_Object &_object, const unsigned &_nb )
        { return lamarck->minimize( _object, _nb);  }
      void print_xmgrace();
      void print_xml()
        { lamarck->print_xml(); };

    protected:
      eoGenOp<t_Object>* make_GenOp(const TiXmlElement *element, std::ofstream &_f);
      eoGenOp<t_Object>* make_recurrent_op(const TiXmlElement &el,
                                           std::ofstream &_f,
                                           std::string &_special,
                                           std::string &_base);
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

