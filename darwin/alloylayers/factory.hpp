//
//  Version: $Id$
//
#ifndef _DARWIN_ALLOY_LAYERS_FACTORY_H_
#define _DARWIN_ALLOY_LAYERS_FACTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/filesystem/path.hpp>
#include <boost/lambda/bind.hpp>
#include <string>
#include <fstream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>
#include <vff/layered.h>

#include "../bitstring.h"
#include "../vff.h"
#include "../bandgap_stubs.h"
#include "../electric_dipole.h"

namespace LaDa
{
  namespace GA
  {
    namespace AlloyLayers
    {
      //! Holds some functions to connect quantities, functionals, printing policies.
      namespace Factory
      {
        //! Connects printing of the genotype with the printing policy.
        template< class T_EVALUATOR > void genotype( const T_EVALUATOR& );
        //! Connects quantities with strain energy, and prints strain energy.
        template< class T_EVALUATOR > void strain_energy( T_EVALUATOR & _evaluator );
        //! Connects epitaxial strain with quantities and printing policies.
        template< class T_EVALUATOR > void epitaxial_strain( T_EVALUATOR &_evaluator );
        //! Connects band-gap with quantities and printing policies
        template< class T_EVALUATOR > void bandgap( T_EVALUATOR &_evaluator );
        //! Connects optical transition strength with quantities and printing policies
        template< class T_EVALUATOR > void transitions( T_EVALUATOR &_evaluator );
        //! Connects valence band offset with quantities and printing policies
        template< class T_EVALUATOR > void vbm( T_EVALUATOR &_evaluator );
        //! Connects conduction band offset with quantities and printing policies
        template< class T_EVALUATOR > void cbm( T_EVALUATOR &_evaluator );
        //! Reads objectives from XML input file.
        template< class T_FACTORY >
          void read_physical_properties( T_FACTORY& _factory, 
                                         boost::filesystem::path &_path );

        //! Prints out connections.
        void declare( const std::string& _string );

        template< class T_EVALUATOR > void genotype( const T_EVALUATOR& )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          t_Object :: connect_print( bl::_1 << bl::ret< LaDa::GA::BitString::Object<> >( bl::_2 )
                                            << bl::constant(" ") ); 
        }

        template< class T_EVALUATOR > void strain_energy( T_EVALUATOR & _evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
          (
            bl::bind
            ( 
              &t_Quantity::push_back,
              bl::_2,
              bl::bind( &t_Object::energy, bl::_1 ) / bl::constant(16.0217733) 
            )
          );
          t_Object :: connect_print
          (
            bl::_1 << bl::constant(" Strain Energy: ")
                   << bl::bind( &t_Object::energy, bl::_2 ) / bl::constant( 16.0217733 )
                   << bl::constant( " eV/f.u. " )
          );
          declare( "Strain energy(VFF)" );
        }

        template< class T_EVALUATOR > void epitaxial_strain( T_EVALUATOR &_evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
          (
            bl::bind
            (
              &t_Quantity::push_back,
              bl::_2,
              bl::bind
              ( 
                &Vff::inplane_stress,
                bl::bind<atat::rMatrix3d>( &t_Object :: stress, bl::_1 ),
                bl::bind<atat::rVector3d>( &T_EVALUATOR :: get_direction,
                                           bl::constant( boost::cref( _evaluator ) ) )
              ) / bl::constant( 16.0217733 )
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
                        bl::bind<atat::rVector3d>( &T_EVALUATOR :: get_direction,
                                                   bl::constant( boost::cref( _evaluator ) ) )
                      ) / bl::constant( 16.0217733 ) << bl::constant( " eV/f.u. " )
          );
          declare( "Epitaxial strain (VFF)" );
        }

        template< class T_EVALUATOR > void bandgap( T_EVALUATOR &_evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
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
            bl::_1 << bl::constant( "Bandgap: ") 
                   << bl::bind( &t_Object::vbm, bl::_2 ) - bl::bind( &t_Object::cbm, bl::_2 )
                   << bl::constant( " " )
          );
          declare( "Band-gaps(Escan)" );
        }

        template< class T_EVALUATOR > void transitions( T_EVALUATOR &_evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.do_dipole( true );
          _evaluator.connect
          (
            bl::bind
            (
              &t_Quantity::push_back,
              bl::_2, bl::bind(  &t_Object::osc_strength, bl::_1 )
            )
          );
          t_Object :: connect_print
          ( 
            bl::_1 << bl::ret<const Keepers::OscStrength&>( bl::_2 )
                   << bl::constant( " " )
          );
          declare( "Transition dipoles(Escan)" );
        }

        template< class T_EVALUATOR > void vbm( T_EVALUATOR &_evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
          (
            bl::bind
            (
              &t_Quantity::push_back,
              bl::_2,
              bl::bind( &t_Object::vbm, bl::_1 )
            )
          );
          t_Object :: connect_print
          ( 
            bl::_1 << bl::constant("VBM: ") << bl::bind( &t_Object::vbm, bl::_2 )
                   << bl::constant( " " )
          );
          declare( "Valence-band offsets (Escan)" );
        }

        template< class T_EVALUATOR > void cbm( T_EVALUATOR &_evaluator )
        {
          namespace bl = boost::lambda;
          typedef typename T_EVALUATOR :: t_Individual :: t_IndivTraits t_IndivTraits;
          typedef typename t_IndivTraits :: t_Object t_Object;
          typedef typename t_IndivTraits :: t_QuantityTraits :: t_Quantity t_Quantity;
          _evaluator.connect
          (
            bl::bind
            (
              &t_Quantity::push_back,
              bl::_2,
              bl::bind( &t_Object::cbm, bl::_1 )
            )
          );
          t_Object :: connect_print
          ( 
            bl::_1 << bl::constant("CBM: ") << bl::bind( &t_Object::cbm, bl::_2 )
                   << bl::constant( " " )
          );
          declare( "Conduction-band offsets (Escan)" );
        }

        void declare( const std::string& _string )
        {
          LaDa::Print::out << _string << " will be optimized.\n";
          LaDa::Print::xmg << LaDa::Print::Xmg::comment
                           <<  _string << " will be optimized.\n";
        }

        template< class T_FACTORY >
          void read_physical_properties( T_FACTORY& _factory,
                                         boost::filesystem::path &_path )
          {
            TiXmlDocument doc; 
            std::string filetxt;
            opt::read_xmlfile( _path, filetxt );
            doc.Parse( filetxt.c_str() );
            TiXmlHandle docHandle( &doc ); 

            const TiXmlElement* child = docHandle.FirstChild( "Job" )
                                                 .FirstChild( "GA" )
                                                 .FirstChild( "Objectives" )
                                                 .FirstChild( "Objective" )
                                                 .Element();
            __DOASSERT( not child, "Could not find Objective tag.\n" )
            for(; child; child = child->NextSiblingElement("Objective") )
            {
              __DOASSERT( not child->Attribute("value"),
                          "Found objective tag without a value attribute.\n" )
              const std::string value = child->Attribute( "value" );
              _factory( value );
            }
          }

      } // namespace factory
    } // namespace Layered

  }
} // namespace LaDa

#include "evaluator.impl.h"

#endif // _LAYERED_H_
