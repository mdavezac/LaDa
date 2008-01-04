//
//  Version: $Id$
//
#ifndef  _DARWIN_H_
#define  _DARWIN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <list>

#include <eo/eoPop.h>
#include <eo/eoGenOp.h>
#include <eo/utils/eoState.h>
#include <eo/eoReplacement.h>

#include <tinyxml/tinyxml.h>

#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

#include "checkpoints.h"
#include "taboos.h"
#include "objective.h"
#include "store.h"
#include "evaluation.h"
#include "breeder.h"
#include "scaling.h"
#include "gatraits.h"

namespace GA
{
  //! \brief Also does your laundry
  //! \details This class is god. It controls the input from XML. It runs the
  //! GA. It ouptuts the results. It calls all the shots. And yet its generic.
  //! There is no space in this class for anything application-specific. Just put
  //! it in \a T_EVALUATOR.
  template< class T_EVALUATOR >
  class Darwin
  {    
    public:
      typedef T_EVALUATOR t_Evaluator; //!< Evaluator type
      //! All %GA traits
      typedef typename t_Evaluator :: t_GATraits            t_GATraits;
    private:
      //! Type of this class
      typedef Darwin<t_Evaluator>                           t_This;
      //! Type of the individuals
      typedef typename t_GATraits :: t_Individual           t_Individual;
      //! %Traits of the quantity (or raw fitness)
      typedef typename t_GATraits :: t_QuantityTraits       t_QuantityTraits;
      //! Type of the genomic object
      typedef typename t_GATraits :: t_Object               t_Object;
      //! Type of the population
      typedef typename t_GATraits :: t_Population           t_Population;
      //! Type of the collection of populations
      typedef typename t_GATraits :: t_Islands              t_Islands;
      //! Type of the scalar quantity (or scalar raw fitness)
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      //! Type of the objective type holder
      typedef typename Objective :: Types < t_GATraits >    t_ObjectiveType;
      //! Type of the storage type holder
      typedef typename Store :: Types< t_GATraits >         t_Store;
      //! Type of the scaling virtual base class
      typedef Scaling :: Base<t_GATraits>                   t_Scaling;

    protected:
      //! \brief Input/Output Flag. \see Darwin::do_save, Darwin::do_restart
      const static types::t_unsigned SAVE_RESULTS     = 1; 
      //! \brief Input/Output Flag. \see Darwin::do_save, Darwin::do_restart
      const static types::t_unsigned SAVE_HISTORY     = 2; 
      //! \brief Input/Output Flag. \see Darwin::do_save, Darwin::do_restart
      const static types::t_unsigned SAVE_POPULATION  = 4; 

    protected:
      //! \brief Filename of input
      //! \details At this point (revision ~310), this is no real option.
      //! Darwin::filename is hardcoded to "input.xml".
      //! \todo allow commandline options for specifying input file
      std::string filename; 
      //! \brief Input file of the evaluator
      //! \details If it is different from Darwin::filename, it should be
      //! indicated in a \<Filename evaluator="?"/\> tag in Darwin::filename. Otherwise,
      //! evaluator input is read from Darwin::filename.
      std::string evaluator_filename;
      //! \brief Input filename where restart data (in XML) can be found
      //! \details It can be different from Darwin::filename. It is set by a
      //! \<Filename restart="?" /\> within the \<GA\> .. \</GA\> tags of
      //! Darwin::filename. Only those components specified by
      //! Darwin::do_restart are read from input.
      std::string restart_filename;
      //! \brief Input filename where restart data (in XML) is saved to
      //! \details It can be different from Darwin::filename. It is set by a
      //! \<Filename restart="?" /\> within the \<GA\> .. \</GA\> tags of
      //! Darwin::filename. Only those components specified by
      //! Darwin::do_save are actually saved to input. 
      //! The following code 
      //! \code
      //  Darwin::do_save & Darwin::SAVE_RESULTS
      //! \endcode
      //! is always true on at least one processor.
      std::string save_filename;
      //! Size of the deterministic tournaments used to choose parents prior to mating
      types::t_unsigned tournament_size;
      //! Size of the population
      types::t_unsigned pop_size;
      //! \brief Maximum of number of generations before quitting.
      //! \details This is not necessarily the only way to quit...
      types::t_unsigned max_generations;
      //! Number of independent islands (eg, independent populations)
      types::t_unsigned nb_islands;
      //! \brief Says which components to save 
      //! \details There must always be at least one processor on which
      //! \code
      //  Darwin::do_save & Darwin::SAVE_RESULTS
      //! \endcode
      //! is always true.
      //! \see Darwin::SAVE_RESULTS, Darwin::SAVE_HISTORY, Darwin::SAVE_POPULATION
      types::t_unsigned do_save;
      //! \brief Says which components to reload
      //! \see Darwin::SAVE_RESULTS, Darwin::SAVE_HISTORY, Darwin::SAVE_POPULATION
      types::t_unsigned do_restart;
      //! ratio of offspring to population size
      types::t_real     replacement_rate;
      //! \brief Print::xmg output flag. 
      //! \see  GA::PrintGA \see Print::xmg
      bool do_print_each_call;
      enum {
        //! Starting population is chosen randomly
        RANDOM_POPULATE,
        //! \brief Starting population is chosen accross a random partition of space
        //! \todo move partition_populate to BitString namespace? 
        //! and make this an application-specific...
        PARTITION_POPULATE 
      } populate_style; //!< Flag for how to choose starting population

      //! \brief The continuator containing all others. See EO library
      //! \see Darwin::Load_CheckPoints()
      IslandsContinuator<t_GATraits>*    continuator;
      //! \brief Pointer to a collection containing previously assessed individuals
      //! \details This pointer can be set to NULL. In that case, no history
      //! tracking is done.
      //! \see Darwin::make_History(), Darwin::Save(), Darwin::Restart()
      History<t_Individual>*             history;
      //! \brief Pointer to the taboo virtual base class
      //! \details This pointer can be set to NULL. In that case, nothing is taboo.
      //! Using a virtual base class allows the instance to be pretty much
      //! anything we want. 
      //! \todo Add restarting and saving capabilities for all taboos in GA::Darwin
      //! \see Darwin::Load_Taboos().
      Taboo_Base<t_Individual>*          taboos;
      //! Mating operators
      //! \see Darwin::make_genetic_op(), Darwin::Load_Mating()
      eoGenOp<t_Individual>*             breeder_ops;
      //! \brief Pointer to a breeder object handling all things mating.
      //! \details This pointer should not be NULL when calling Darwin::run. 
      //! \todo make Darwin::breeder an instance rather than an object?
      //! \see Darwin::Load_Mating(), Darwin::make_breeder()
      Breeder<t_GATraits>*               breeder;
      //! \brief Pointer to the replacement scheme
      //! \details The replacement shemes say how to go from an old population
      //! and a (previously created) offspring population to a new population. Eg
      //! it does culling and meshing operations. This pointer should not be
      //! NULL when calling Darwin::run(). 
      //! \see Darwin::make_breeder()
      eoReplacement<t_Individual>*       replacement;
      //! \brief Pointer to an objective
      //! \details This objective is vectorial by default in  multi-objective
      //! applications. It is scalar in single-objective applications. 
      //! \see Darwin::Load_Method()
      typename t_ObjectiveType::Vector*  objective;
      //! \brief Pointer to a class implementing storage capabilities
      //! \see Darwin::Load_Method()
      typename t_Store :: Base*          store;
      //! \brief Pointer to a funnel through which all evaluations are done
      //! \see Darwin::Load_Method()
      Evaluation::Base<t_GATraits>*      evaluation;
      //! \brief Points to an object which applies population-dependent
      //! operations on individual fitnesses
      //! \details  This pointer can be set to NULL. In that case, nothing is taboo.
      //! Using a virtual base class allows the instance to be pretty much
      //! anything we want. 
      //! \see Darwin::Load_Method()
      t_Scaling*                         scaling;

      //! \brief The evaluator instance itself
      //! \see Darwin::Load(const TiXmlElement& )
      t_Evaluator   evaluator;
      //! holds a collection of independent populations
      t_Islands     islands;
      //! Offspring population
      t_Population  offspring;
      
      //! \brief Contains most, if not all, pointers allocated in Darwin. 
      //! \details Is an EO benefit. Basically takes care of deallocating all
      //! those pointers once the game is up. Go check out EO to find out how it
      //! works.
      eoState eostates;

    public:
      //! Constructor
      Darwin () : filename("input.xml"), tournament_size(2), pop_size(100),
                  max_generations(0), nb_islands(1), do_save(SAVE_RESULTS),
                  do_restart(0), replacement_rate(0.1), do_print_each_call(false),
                  populate_style(RANDOM_POPULATE), continuator(NULL), history(NULL),
                  taboos(NULL), breeder_ops(NULL), breeder(NULL), replacement(NULL),
                  objective(NULL), store(NULL), evaluation(NULL), scaling(NULL) {}
      //! Destructor
      virtual ~Darwin ();

      //! \brief Load the XML input file \a _filename
      //! \details Directs all the loading for the genetic algorithm. It makes
      //! all the calls, from %GA input to requesting a load from the evaluator,
      //! to restarting history... Basically, if you can't find it somewhere
      //! from here, it's probably not gonna get loaded.
      //! \note This %%function makes a call to Evaluator::Load(const
      //!       TiXmlElement&) and Evaluator::Load(std::string&)
      bool Load( const std::string &_filename );
      //! \brief Runs the generational loop
      //! \brief Executive member of the genetic algorithm. It starts by
      //! initializing a population, and goes on to evolving through the
      //! generational loop. Much as Darwin::Load, if you can't find it called
      //! somewhere from here, it's probly not gonna get exectuted.
      //! Nonetheless, the loop itself is rather simple, so go look at it.
      void run();

    protected: 
      //! \brief Loads overall Genetic Algorithm attributes
      //! \note This %function makes a call to Evaluator::LoadAttribute()
      bool Load_Parameters( const TiXmlElement &_parent );
      //! Creates the history object, if requested on input
      void make_History( const TiXmlElement &_parent );
      //! Creates the mating operations from input
      bool Load_Mating( const TiXmlElement &_parent );
      //! \brief Creates the objectives, the storage, the evaluation, all that says
      //! what it is you are looking for. 
      void Load_Method( const TiXmlElement &_parent );
      //! \brief Creates taboo objects if required on input
      //! \note This %function makes a call to Evaluator::LoadTaboo()
      void Load_Taboos( const TiXmlElement &_node );
      //! Creates checkpoints, such as for printing, saving...
      //! \note This %function makes a call to Evaluator::LoadContinue()
      void Load_CheckPoints (const TiXmlElement &_parent);
      //! Loads "Restart" information from an XML node \a _node
      bool Restart(const TiXmlElement &_node);
      //! Loads "Restart" information from file
      bool Restart();
      //! Saves "Restart" information to an XML node \a _node
      bool Save(const TiXmlElement &_node);
      //! Saves "Restart" information to file
      bool Save();
      //! \brief Creates genetic operators from input
      //! \return a pointer to the genetic operator specified by \a el. It
      //!         returns NULL, if the XML input is incorrect, or not a genetic
      //!         operator.
      //! \details This %function is capable of recursively calling itself. In
      //!          other words it will automatically create a genetic operator
      //!          containing other genetic operators. In fact, it will keep
      //!          going 'till it has exhausted the meaning of \a el. Just
      //!          call it without specifying the \a current_op argument, eg
      //!          call the one argument version and store the return value.
      //! \note This %function makes a call to Evaluator::LoadGaOp()
      //! \param el XML node at the start of a genetic operator
      //! \param current_op Don't touch, don't use, and ignore. This is for
      //!                   internal use only.
      eoGenOp<t_Individual>* make_genetic_op( const TiXmlElement &el,
                                              eoGenOp<t_Individual> *current_op = NULL);
      //! Make the breeder object
      void make_breeder();
      //! Creates the replacement scheme
      eoReplacement<t_Individual>* make_replacement();

      //! \brief Creates a new starting population
      //! \see Darwin::random_populate(), Darwin::partition_populate()
      void populate ();
      //! \brief Creates a new random population
      //! \note This %function makes a call to Evaluator::initialize()
      void random_populate ( t_Population &_pop, types::t_unsigned _size);
      //! \brief Creates a new random population using a partition scheme
      //! \note This %function makes a call to Evaluator::initialize() and Object::mask()
      void partition_populate ( t_Population &_pop, types::t_unsigned _size);

      //! \brief Submits individuals to history, taboo, etc, prior to starting %GA
      //! \details initializes the endopoints of a convex-hull, for instance.
      //! Presubmitted individuals are not put inot the population.
      //! \note This %function makes a call to Evaluator::presubmit()
      void presubmit();
#ifdef _MPI
      /** \ingroup MPI
       * \brief Loads the input filenames and input files, and broadcasts them
       * to all processors.
       * \details  It is quite a pain to write broadcasting functions which
       * will correctly create the mostly virtual and imbricated types we have
       * in %GA. Its much easier to have each processor process the input. By
       * reading the input files as a whole and broadcasting the results, we
       * can process the input in parallel. 
       */
      void LoadAllInputFiles(std::string &_input, 
                             std::string &_restart, 
                             std::string &_evaluator );
#endif
  };

// template< class T_EVALUATOR >
//   const types::t_unsigned Darwin<T_EVALUATOR> :: SAVE_RESULTS    = 1;
// template< class T_EVALUATOR >
//   const types::t_unsigned Darwin<T_EVALUATOR> :: SAVE_POPULATION = 2;
// template< class T_EVALUATOR >
//   const types::t_unsigned Darwin<T_EVALUATOR> :: SAVE_HISTORY     = 4;

} // namespace GA

#include "darwin.impl.h"

#endif // _DARWIN_H_

