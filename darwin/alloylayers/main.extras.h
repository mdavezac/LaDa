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
#   define __PROGNAME__ "alloylayers_opt"

    typedef Individual :: Types
            < 
              GA :: AlloyLayers :: Object,
              TwoSites :: Concentration,
              TwoSites :: Fourier 
            > :: t_Vector t_Individual;
    typedef GA :: AlloyLayers :: Evaluator
            <
              t_Individual,
              GA :: AlloyLayers :: Translate,
              GA :: AlloyLayers :: AssignSignal
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
         "Specify properties to optimize on input. Call --help.\n" 
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
    Print::out << "Strain energy (Vff) will"
               << ( optimize_strain ? " ": " NOT " )
               << " be optimized.\n";
    Print::xmg << Print::Xmg::comment
               << "Strain energy (Vff) will"
               << ( optimize_strain ? " ": " NOT " )
               << " be optimized." << Print::endl;
    Print::out << "Epitaxial strain (Vff) will"
               << ( optimize_epi ? " ": " NOT " )
               << " be optimized.\n";
    Print::xmg << Print::Xmg::comment
               << "Epitaxial strain (Vff) will"
               << ( optimize_epi ? " ": " NOT " )
               << " be optimized." << Print::endl;
    Print::out << "Band-gaps (escan) will"
               << ( optimize_bandgaps ? " ": " NOT " )
               << " be optimized.\n";
    Print::xmg << Print::Xmg::comment
               << "Band-gaps (escan) will"
               << ( optimize_bandgaps ? " ": " NOT " )
               << " be optimized." << Print::endl;
    Print::out << "Optical transitions (escan/oscillator strength) will"
               << ( optimize_transitions ? " ": " NOT " )
               << " be optimized.\n";
    Print::xmg << Print::Xmg::comment
               << "Optical transitions (escan/oscillator strength) will"
               << ( optimize_transitions ? " ": " NOT " )
               << " be optimized." << Print::endl;

# elif _MAIN_ALLOY_LAYERS_EXTRAS_ == 3
    // Connects assignement functors and print functors depending on requested properties.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_ 4

    typedef t_Evaluator :: t_GATraits :: t_QuantityTraits :: t_Quantity t_Quantity;
    typedef t_Evaluator :: t_GATraits :: t_Object t_Object;
    // Prints object.
    if( print_genotype )
      t_Object :: connect_print( bl::_1 << bl::ret< BitString::Object<> >( bl::_2 ) );
    if( optimize_strain )
    {
      ga.evaluator.connect
      (
        bl::bind( &t_Quantity::push_back,
                  bl::_2, bl::bind( &t_Object::energy, bl::_1 ) )
      );
      t_Object :: connect_print
      (
        bl::_1 << bl::constant(" Strain Energy: ")
               << bl::bind( &t_Object::energy, bl::_2 )
      );
    }
    if( optimize_epi )
    {
      ga.evaluator.connect
      (
        bl::bind
        (
          &t_Quantity::push_back,
          bl::_2,
          bl::bind
          ( 
            &Vff::inplane_stress,
            bl::bind<atat::rMatrix3d>( &t_Object :: stress, bl::_1 ),
            bl::constant( ga.evaluator.get_direction() )
          )
        )
      );
      t_Object :: connect_print
      (
        bl::_1 << bl::constant(" Epi Strain: ")
               << bl::constant( std::fixed ) << bl::constant( std::setw(8) )
               << bl::constant( std::setprecision(3) )
               << bl::bind
                  ( 
                    &Vff::inplane_stress,
                    bl::bind<atat::rMatrix3d>( &t_Object :: stress, bl::_2 ),
                    bl::constant( ga.evaluator.get_direction() )
                  ) / 16.0217733 << bl::constant( " eV/structure" )
      );
    }
    if( optimize_bandgaps )
    {
      ga.evaluator.connect
      (
        bl::bind
        (
          &t_Quantity::push_back,
         bl::_2,
          bl::bind( &t_Object::vbm, bl::_1 ) - bl::bind( &t_Object::cbm, bl::_1 )
        )
      );
      t_Object :: connect_print
      ( 
        bl::_1 << bl::ret<const ::GA::Keepers::BandGap&>( bl::_2 )
      );
    }
    ga.evaluator.do_dipole = optimize_transitions;
    if( optimize_transitions )
    {
      ga.evaluator.connect
      (
        bl::bind
        (
          &t_Quantity::push_back,
          bl::_2, bl::bind(  &t_Object::osc_strength, bl::_1 )
        )
      );
      t_Object :: connect_print
      ( 
        bl::_1 << bl::ret<const ::GA::Keepers::OscStrength&>( bl::_2 )
      );
    }

# else 
    // makes sure this file cannot be (meaningfully) included anymore.
#   undef _MAIN_ALLOY_LAYERS_EXTRAS_
#   define _MAIN_ALLOY_LAYERS_EXTRAS_
# endif

#endif
