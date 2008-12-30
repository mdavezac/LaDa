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
#   include <factory/factory.h>
#   include "factory.hpp"
#   include "../operators/periodic.h"
#   include "../operators/factory.h"
#   include "../operators/xmlfactory.h"
#   include "../operators/eogenop_adapter.h"
#   include "../operators/taboo.h"
#   include "../bitstring/mutation.h"
#   include "../bitstring/crossover.h"
#   include "operators.h"
#   include "../assign_callbacks.h"
#   include "../static_translate.h"
#   define __PROGNAME__ "Epitaxial Alloy Layer Optimization"

    typedef LaDa::Individual :: Types
            < 
              LaDa::GA::static_translate
              <
                LaDa :: GA :: AlloyLayers :: Object,
                LaDa :: GA :: AlloyLayers :: Translate
              > :: type,
              LaDa :: TwoSites :: Concentration,
              LaDa :: TwoSites :: Fourier 
            > :: t_Vector t_Individual;
    typedef LaDa :: GA :: AlloyLayers :: Evaluator
            <
              t_Individual,
              LaDa :: GA :: AlloyLayers :: Translate,
              LaDa :: GA :: AssignSignal
            > t_Evaluator;

   template< class T_TYPE > T_TYPE call_lexcast( const std::string& _str )
    { return boost::lexical_cast<T_TYPE>( _str ); }

   inline void reference_updater( bool, LaDa::BandGap::Darwin<LaDa::Vff::Layered> & _bg )
     { _bg.Continue(); }

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
    ga.evaluator.do_dipole( false );
    LaDa::Factory::Factory<void(void), std::string> properties_factory;
    properties_factory.connect
      ( "epi", "Epitaxial Strain  in eV per f.u.\n"
        "(trace of strain - strain along growth direction).",
        boost::bind<void>( &PPFactory::epitaxial_strain<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "energy", "VFF energy in eV per f.u." ,
        boost::bind<void>( &PPFactory::strain_energy<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "bandgap", "Band gap in eV", 
        boost::bind( &PPFactory::bandgap<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "transition", "Dipole oscillator strength between VBM and CBM. Arbitrary units.",
        boost::bind( &PPFactory::transitions<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "cbm", "Conduction band minimum in eV.", 
        boost::bind( &PPFactory::cbm<t_Evaluator>, boost::ref(ga.evaluator) ) )
      ( "vbm", "Valence band minimum in eV.",
        boost::bind( &PPFactory::vbm<t_Evaluator>, boost::ref(ga.evaluator) ) );

    // connects ga operator factory.
    ga.operator_factory.connect
      ( "Random", "creates a random individual.",
        boost::bind( &PPFactory::random<t_OpFactory, t_Evaluator>,
                     _1, _2, _3, boost::ref(ga.evaluator) ) )
      ( "Crossover", "bitstring-like crossover.",
        boost::bind( &LaDa::GA::BitString::crossover_factory<t_OpFactory>, _1, _2, _3 ) )
      ( "Mutation", "bistring-like mutation.",
        boost::bind( &LaDa::GA::BitString::mutation_factory<t_OpFactory>, _1, _2, _3 ) );


# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 4
    // reads program options.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 5
    if( print_genotype ) PPFactory::genotype( ga.evaluator );
    PPFactory::read_physical_properties( properties_factory, input );
    ga.checkpoints.connect_updater( boost::bind( &reference_updater, _1, 
                                                 boost::ref( ga.evaluator.bandgap ) ) );

# else 
    // makes sure this file cannot be (meaningfully) included anymore.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ nofurtherinclude
# endif

#endif
