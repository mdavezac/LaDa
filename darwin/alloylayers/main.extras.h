//
//  Version: $Id$
//


// This file contains random code snipets for _ALLOY_LAYERS_
#ifdef _ALLOY_LAYERS_

// includes.
# if ! defined( _MAIN_ALLOY_LAYERS_EXTRAS_ )
    // includes some more files 
    // defines individual and evaluator.
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 0
#   include <eo/eoPopulator.h>
#   include <boost/lambda/bind.hpp>
#   include <boost/bind.hpp>
#   include "../vff.h"
#   include "../two_sites.h"
#   include "evaluator.h"
#   include "policies.h"
#   include "object.h"
#   include <opt/factory.h>
#   include "factory.hpp"
#   include "../operators/periodic.h"
#   include "../operators/factory.h"
#   include "../operators/xmlfactory.h"
#   include "../operators/eogenop_adapter.h"
#   include "operators.h"
#   include "../taboos.h"
#   define __PROGNAME__ "Epitaxial Alloy Layer Optimization"

    typedef LaDa::Individual :: Types
            < 
              LaDa :: GA :: AlloyLayers :: Object,
              LaDa :: TwoSites :: Concentration,
              LaDa :: TwoSites :: Fourier 
            > :: t_Vector t_Individual;
    typedef LaDa :: GA :: AlloyLayers :: Evaluator
            <
              t_Individual,
              LaDa :: GA :: AlloyLayers :: Translate,
              LaDa :: GA :: AlloyLayers :: AssignSignal
            > t_Evaluator;

   template< class T_TYPE > T_TYPE call_lexcast( const std::string& _str )
    { return boost::lexical_cast<T_TYPE>( _str ); }

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 0
    // defines program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 1

      ("printg,g", "Print genotype when printing invidivual." )

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 1
    // reads program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 2

      const bool print_genotype = vm.count("printg") > 0;
   //  std::cout << "Bitstring genotype will"
   //            << ( print_genotype ? " ": " NOT " )
   //            << " be printed.\n";

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 2 || _MAIN_ALLOY_LAYERS_EXTRAS_ == 3
    // Creates factories.
#   if _MAIN_ALLOY_LAYERS_EXTRAS_ == 2
#     undef _MAIN_ALLOY_LAYERS_EXTRAS_
#     define _MAIN_ALLOY_LAYERS_EXTRAS_ 3
#   else
#     undef _MAIN_ALLOY_LAYERS_EXTRAS_
#     define _MAIN_ALLOY_LAYERS_EXTRAS_ 4
#   endif

    namespace bl = boost::lambda;
    short_description = "Performs optimizations within the  \"alloy-layers\""
                        " configuration space of Zinc-Blend structures.\n";
    namespace PPFactory = LaDa::GA::AlloyLayers::Factory;
    typedef LaDa :: GA :: Darwin< t_Evaluator > t_Darwin;

    // connects physical properties factory.
//   ga.evaluator.do_dipole( false );
    LaDa::Factory::Factory<void(void), std::string> properties_factory;
    properties_factory.connect
      ( "epi", "Epitaxial Strain  in eV per f.u.\n"
        "          (trace of strain - strain along growth direction).",
        boost::bind<void>( &PPFactory::epitaxial_strain<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "energy", "VFF energy in eV per f.u." ,
        boost::bind<void>( &PPFactory::strain_energy<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "bandgap", boost::bind( &PPFactory::bandgap<t_Evaluator>, boost::ref(ga.evaluator) ) )
//     ( "transition", "Dipole oscillator strength between VBM and CBM. Arbitrary units.",
//       boost::bind( &PPFactory::transitions<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "cbm", "Conduction band minimum in eV.", 
        boost::bind( &PPFactory::cbm<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "vbm", "Valence band minimum in eV.",
        boost::bind( &PPFactory::vbm<t_Evaluator>, boost::ref(ga.evaluator) ) );

    // connects ga operator factory.
    typedef LaDa :: GA :: Factory :: XmlOperators< t_Individual, eoPopulator<t_Individual> >
      t_OpFactory;
    t_OpFactory op_factory;
    op_factory.connect_attribute
      ( "period", "operators with this attribute will be\n"
        "             called only every n=\"period\" generations.\n",
         boost::bind( &LaDa::GA::Factory::periodic< t_OpFactory >,
                       _1, _2, _3, ga.get_counter() ) );
    op_factory.connect
      ( "Operators", "container of operators. Accepts attribute "
                     " type=\"and\" or type=\"or\".",
        boost::bind( &LaDa::GA::Factory::containers<t_OpFactory>, _1, _2, _3 ) ) 
      ( "And", "Sequential container. same as <Operators type=\"and\"> ... <Operators/>.",
        boost::bind( &LaDa::GA::Factory::sequential<t_OpFactory>, _1, _2, _3 ) ) 
      ( "Or", "Proportional container. same as <Operators type=\"or\"> ... <Operators/>.",
        boost::bind( &LaDa::GA::Factory::sequential<t_OpFactory>, _1, _2, _3 ) )
      ( "TabooOp", "An offspring created by the operators within this one\n"
        "              will be rejected if it is in the tabooed individuals \n"
        "              defined in <Taboos> .. </Taboos>. Offspring creation\n"
        "              will repeat for at most n (attribute max=\"n\") trials,\n"
        "              at which point a random individual is created.",
        boost::bind( &LaDa::GA::Factory::taboo_op<t_OpFactory, t_Darwin>, 
                     _1, _2, _3, boost::ref(ga), "Random" ) )
      ( "Random", "creates a random individual.",
        boost::bind( &PPFactory::random<t_OpFactory, t_Evaluator>,
                     _1, _2, _3, boost::ref(ga.evaluator) ) )
      ( "Crossover", "bitstring-like crossover.",
        boost::bind( &PPFactory::crossover<t_OpFactory>, _1, _2, _3 ) )
      ( "Mutation", "bistring-like mutation.",
        boost::bind( &PPFactory::mutation<t_OpFactory>, _1, _2, _3 ) );
    ga.set_operator_factory( op_factory );
    if( print_genotype ) PPFactory::genotype( ga.evaluator );
    PPFactory::read_physical_properties( properties_factory, input );

    LaDa::Factory::Factory< void(const std::string& ), std::string >& att_factory = ga.att_factory;
    ga.tournament_size = 2;
    ga.replacement_rate = 0.5;
    ga.pop_size = 1;
    ga.max_generations = 0;
    att_factory.connect
      (
        "tournament", "size of the deterministic tournaments\n"
        "           for each parent selection. Default = 2.",
        bl::var(ga.tournament_size) = bl::bind( &call_lexcast<LaDa::types::t_unsigned>, bl::_1 )
      )
      (
        "rate", "      replacement rate (eg offspring creation rate).\n"
        "           Default = 0.5.",
        bl::var(ga.replacement_rate) = bl::bind( &call_lexcast<LaDa::types::t_real>, bl::_1 ) 
      )
      (
        "popsize", "   population size. Default = 1.",
        bl::var(ga.pop_size) = bl::bind( &call_lexcast<LaDa::types::t_unsigned>, bl::_1 ) 
      )
      (
        "maxgen", "    Max number of generations before quitting."
        "           Default = 0.",
        bl::var(ga.max_generations) = bl::bind( &call_lexcast<LaDa::types::t_unsigned>, bl::_1 ) 
      )
      (
        "islands", "    Number of independent islands. Default = 1.",
        bl::var(ga.nb_islands) = bl::bind( &call_lexcast<LaDa::types::t_unsigned>, bl::_1 ) 
      )
      (
        "seed", "       random seed. Default = 0."
        "            In MPI, seed0, seed1, ... are accepted,"
        "            if not necessarily used.",
        boost::bind( &LaDa::GA::Topology::seed_n, ga.topology, 0, _1 )
      );
#    ifdef _MPI
        for( size_t ii(0); ii < LaDa::mpi::main->size(); ++ii )
          att_factory.connect
          (
            "seed" + boost::lexical_cast<std::string>(ii),
            "       random seed. Default = 0.",
            boost::bind( &LaDa::GA::Topology::seed_n, ga.topology, 0, _1 )
          );
#     endif

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 4
    // reads program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 5
# else 
    // makes sure this file cannot be (meaningfully) included anymore.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_
# endif

#endif
