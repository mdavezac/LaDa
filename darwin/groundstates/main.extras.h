//
//  Version: $Id$
//


// This file contains random code snipets for Cluster Expansion Optimization.
#ifdef _CE

// includes.
# if ! defined( _MAIN_CE_EXTRAS_ )
    // includes some more files 
    // defines individual and evaluator.
#   define _MAIN_CE_EXTRAS_ 0
#   include <eo/eoPopulator.h>
#   include <boost/lambda/bind.hpp>
#   include <boost/bind.hpp>
#   include "../single_site.h"
#   include "evaluator.h"
#   include "object.h"
#   include <factory/factory.h>
#   include "../operators/periodic.h"
#   include "../operators/factory.h"
#   include "../operators/xmlfactory.h"
#   include "../operators/taboo.h"
#   include "../bitstring/mutation.h"
#   include "../bitstring/random.h"
#   include "../bitstring/crossover.h"
#   include "../static_translate.h"
#   include "krossover.h"
#   include "../single_site.h"
#   include "populate_factory.h"
#   define __PROGNAME__ "Cluster Expansion Optimization"

    typedef LaDa::Individual :: Types
            < 
              LaDa :: GA :: static_translate 
              <
                LaDa :: GA :: GroundStates :: Object,
                LaDa :: GA :: GroundStates :: Translate
              > :: type,
              LaDa :: SingleSite :: Concentration,
              LaDa :: SingleSite :: Fourier 
            > :: t_Scalar t_Individual;
    typedef LaDa :: GA :: GroundStates :: Evaluator
            <
              t_Individual,
              LaDa :: GA :: GroundStates :: Translate,
              LaDa :: GA :: GroundStates :: AssignCE
            > t_Evaluator;

   template< class T_TYPE > T_TYPE call_lexcast( const std::string& _str )
    { return boost::lexical_cast<T_TYPE>( _str ); }

# elif _MAIN_CE_EXTRAS_ == 0
    // defines program options.
#   undef _MAIN_CE_EXTRAS_
#   define _MAIN_CE_EXTRAS_ 1

      ("printg,g", "Print genotype when printing invidivual." )

# elif _MAIN_CE_EXTRAS_ == 1
    // reads program options.
#   undef _MAIN_CE_EXTRAS_
#   define _MAIN_CE_EXTRAS_ 2

      const bool print_genotype = vm.count("printg") > 0;
   //  std::cout << "Bitstring genotype will"
   //            << ( print_genotype ? " ": " NOT " )
   //            << " be printed.\n";

# elif _MAIN_CE_EXTRAS_ == 2 || _MAIN_CE_EXTRAS_ == 3
    // Creates factories.
#   if _MAIN_CE_EXTRAS_ == 2
#     undef _MAIN_CE_EXTRAS_
#     define _MAIN_CE_EXTRAS_ 3
#   else
#     undef _MAIN_CE_EXTRAS_
#     define _MAIN_CE_EXTRAS_ 4
#   endif

    namespace bl = boost::lambda;
    short_description = "Performs optimizations for cluster-expansion functionals.\n";
    typedef LaDa :: GA :: Darwin< t_Evaluator > t_Darwin;

    // connects ga operator factory.
    ga.operator_factory.connect
      ( "Random", "creates a random individual.",
        boost::bind( &LaDa::GA::BitString::random_two_valued_factory<t_OpFactory>,
                     _1, _2, _3, 1, 0 ) )
      ( "Krossover", "k-space crossover.",
        boost::bind( &LaDa::GA::GroundStates::krossover_factory<t_OpFactory>, _1, _2, _3,
                     boost::ref( ga.evaluator.get_structure()) ) )
      ( "Crossover", "bitstring crossover.",
        boost::bind( &LaDa::GA::BitString::crossover_factory<t_OpFactory>, _1, _2, _3 ) )
      ( "Mutation", "bistring mutation.",
        boost::bind( &LaDa::GA::BitString::mutation_factory<t_OpFactory>, _1, _2, _3 ) );
    ga.att_factory.connect
      (
        "populate", "To creation style for the initial population: \n"
        "   = \"random\" (default),\n"
        "   = \"partition\" a few random, plus their photo-negatives,"
        "and so on.",
        boost::bind
        (
          &LaDa::GA::GroundStates::populate_factory<t_Evaluator>,
          _1, boost::ref(ga.population_creator),
          boost::ref(ga.evaluator), boost::cref(ga.taboos) 
        )
      );

# elif _MAIN_CE_EXTRAS_ == 4
    // reads program options.
#   undef _MAIN_CE_EXTRAS_
#   define _MAIN_CE_EXTRAS_ 5

# else 
    // makes sure this file cannot be (meaningfully) included anymore.
#   undef _MAIN_CE_EXTRAS_
#   define _MAIN_CE_EXTRAS_ nofurtherinclude
# endif

#endif
