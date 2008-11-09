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

    namespace factory = LaDa::GA::AlloyLayers::Factory;
    ga.evaluator.do_dipole = false;
    if( print_genotype ) factory::genotype( ga.evaluator );
    LaDa::Factory::PureCalls properties_factory;
    properties_factory.connect
      ( "strain", boost::bind( &factory::strain_energy<t_Evaluator>, ga.evaluator ) )
      ( "epi", boost::bind( &factory::epitaxial_strain<t_Evaluator>, ga.evaluator ) )
      ( "bandgap", boost::bind( &factory::bandgap<t_Evaluator>, ga.evaluator ) )
      ( "transition", boost::bind( &factory::transitions<t_Evaluator>, ga.evaluator ) )
      ( "cbm", boost::bind( &factory::cbm<t_Evaluator>, ga.evaluator ) )
      ( "vbm", boost::bind( &factory::vbm<t_Evaluator>, ga.evaluator ) );
    factory::read_physical_properties( properties_factory, input );

# else 
    // makes sure this file cannot be (meaningfully) included anymore.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_
# endif

#endif
