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

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 0
    // defines program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 1

      ("printg,g", "Print genotype when printing invidivual." )
      ("properties,p", "Lists physical properties available "
                       "for optimization and exits." )

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 1
    // reads program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 2

      const bool print_properties = vm.count("properties") > 0;
      if( print_properties )
      {
        std::cout <<
          "Physical properties available for optimization:\n"
          "  _ Strain Energy           =>   <Objective value=\"strain\"/>\n"
          "  _ Epitaxial Strain        =>   <Objective value=\"epi\"/>\n"
          "  _ Bandgap                 =>   <Objective value=\"bandgap\"/>\n"
          "  _ Valence-band offset     =>   <Objective value=\"vbm\"/>\n"
          "  _ Conduction-band offset  =>   <Objective value=\"cbm\"/>\n"
          "  _ Oscillator strength     =>   <Objective value=\"transition\"/>\n\n";
        return 0;
      }
      const bool print_genotype = vm.count("printg") > 0;
      std::cout << "Bitstring genotype will"
                << ( print_genotype ? " ": " NOT " )
                << " be printed.\n";

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 2
    // Connects assignement functors and print functors depending on requested
    // properties.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 3

    namespace PPFactory = LaDa::GA::AlloyLayers::Factory;
    ga.evaluator.do_dipole = false;
    if( print_genotype ) PPFactory::genotype( ga.evaluator );
    LaDa::Factory::PureCalls properties_factory;
    properties_factory.connect
      ( "strain", boost::bind( &PPFactory::strain_energy<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "epi", boost::bind( &PPFactory::epitaxial_strain<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "bandgap", boost::bind( &PPFactory::bandgap<t_Evaluator>,
                                boost::ref(ga.evaluator) ) )
      ( "transition", boost::bind( &PPFactory::transitions<t_Evaluator>,
                                   boost::ref(ga.evaluator) ) )
      ( "cbm", boost::bind( &PPFactory::cbm<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "vbm", boost::bind( &PPFactory::vbm<t_Evaluator>, boost::ref(ga.evaluator) ) );
    typedef LaDa :: GA :: Factory :: XmlOperators< t_Individual > t_OpFactory;
    t_OpFactory op_factory;
    op_factory.connect_attribute
      ( "period", boost::bind( &LaDa::GA::Factory::periodic< t_OpFactory >,
                               _1, _2, _3, ga.get_counter() ) );
    op_factory.connect
      ( "Operators", boost::bind( &LaDa::GA::Factory::containers<t_OpFactory>, _1, _2, _3 ) )
      ( "Sequential", boost::bind( &LaDa::GA::Factory::sequential<t_OpFactory>, _1, _2, _3 ) )
      ( "Proportional", boost::bind( &LaDa::GA::Factory::sequential<t_OpFactory>, _1, _2, _3 ) )
      ( "TabooOp", boost::bind( &LaDa::GA::Factory::taboo_op<t_OpFactory>, _1, _2, _3,
                                boost::bind( &LaDa::GA::Darwin<t_Evaluator>::get_taboos,
                                             boost::ref( ga ) ), "Random" ) )
      ( "Random", boost::bind( &PPFactory::random<t_OpFactory, t_Evaluator>, _1, _2, _3,
                               boost::ref(ga.evaluator) ) )
      ( "Crossover", boost::bind( &PPFactory::crossover<t_OpFactory>, _1, _2, _3 ) )
      ( "Mutation", boost::bind( &PPFactory::mutation<t_OpFactory>, _1, _2, _3 ) );


# else 
    // makes sure this file cannot be (meaningfully) included anymore.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_
# endif

#endif
