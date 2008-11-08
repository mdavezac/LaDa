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

    ("strain,s", "Optimize strain energy" )
    ("epi,e", "Optimize epitaxial strain" )
    ("bandgaps,b", "Optimize bandgaps" )
    ("transitions,t", "Optimize transition moments" )
    ("printg,g", "Print genotype when printing invidivual." )

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 1
    // reads program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 2


      const bool optimize_strain = vm.count("strain") > 0;
      const bool optimize_epi = vm.count("epi") > 0;
      const bool optimize_bandgaps = vm.count("bandgaps") > 0;
      const bool optimize_transitions = vm.count("transitions") > 0;
      const bool print_genotype = vm.count("printg") > 0;
      __DOASSERT
      ( 
         not ( optimize_strain or optimize_epi or optimize_bandgaps or optimize_transitions ),
         "Specify properties to optimize on input.\nSee " << argv[0] << " --help.\n\n" 
      )
      std::cout << "Strain energy (Vff) will"
                << ( optimize_strain ? " ": " NOT " )
                << " be optimized.\n";
      std::cout <<  "Epitaxial strain (Vff) will"
                <<  ( optimize_epi ? " ": " NOT " )
                <<  " be optimized.\n";
      std::cout << "Band-gaps (escan) will"
                << ( optimize_bandgaps ? " ": " NOT " )
                << " be optimized.\n";
      std::cout << "Optical transitions (escan/oscillator strength) will"
                << ( optimize_transitions ? " ": " NOT " )
                << " be optimized.\n";
      std::cout << "Bitstring genotype will"
                << ( print_genotype ? " ": " NOT " )
                << " be printed.\n";

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 2
    // reads program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 3
      LaDa::Print::out << "Epitaxial strain (Vff) will"
                       << ( optimize_epi ? " ": " NOT " )
                       << " be optimized.\n";
      LaDa::Print::xmg << LaDa::Print::Xmg::comment
                       << "Epitaxial strain (Vff) will"
                       << ( optimize_epi ? " ": " NOT " )
                       << " be optimized." << LaDa::Print::endl;
      LaDa::Print::out << "Band-gaps (escan) will"
                       << ( optimize_bandgaps ? " ": " NOT " )
                       << " be optimized.\n";
      LaDa::Print::xmg << LaDa::Print::Xmg::comment
                       << "Band-gaps (escan) will"
                       << ( optimize_bandgaps ? " ": " NOT " )
                       << " be optimized." << LaDa::Print::endl;
      LaDa::Print::out << "Optical transitions (escan/oscillator strength) will"
                       << ( optimize_transitions ? " ": " NOT " )
                       << " be optimized.\n";
      LaDa::Print::xmg << LaDa::Print::Xmg::comment
                       << "Optical transitions (escan/oscillator strength) will"
                       << ( optimize_transitions ? " ": " NOT " )
                       << " be optimized." << LaDa::Print::endl;

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 3
    // Connects assignement functors and print functors depending on requested properties.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 4

    namespace factory = LaDa::GA::AlloyLayers::Factory;
    ga.evaluator.do_dipole = false;
    if( print_genotype ) factory::genotype( ga.evaluator );
    LaDa::Factory::PureCalls objectives_factory;
    objectives_factory.connect
      ( "strain", boost::bind( &factory::strain_energy<t_Evaluator>, ga.evaluator ) )
      ( "epi", boost::bind( &factory::epitaxial_strain<t_Evaluator>, ga.evaluator ) )
      ( "bandgap", boost::bind( &factory::bandgap<t_Evaluator>, ga.evaluator ) )
      ( "transition", boost::bind( &factory::transition<t_Evaluator>, ga.evaluator ) )
      ( "cbm", boost::bind( &factory::cbm<t_Evaluator>, ga.evaluator ) )
      ( "vbm", boost::bind( &factory::vbm<t_Evaluator>, ga.evaluator ) )
    factory::read_objectives( objectives_factory, input );

# else 
    // makes sure this file cannot be (meaningfully) included anymore.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_
# endif

#endif
