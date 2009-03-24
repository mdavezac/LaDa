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
#include <iomanip>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>
#include <vff/layered.h>
#include <print/stdout.h>
#include <print/xmg.h>

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
        namespace details
        {
          //! Actually function which does the connecting.
          template< class T_EVALUATOR, class T_PROPERTY >  class ConnectProperty;
        }
        //! Connects printing of the genotype with the printing policy.
        template< class T_EVALUATOR > void genotype( const T_EVALUATOR& );
        //! Connects quantities with strain energy, and prints strain energy.
        template< class T_EVALUATOR > void strain_energy( T_EVALUATOR & _evaluator );
        //! Connects epitaxial strain with quantities and printing policies.
        template< class T_EVALUATOR > void epitaxial_strain( T_EVALUATOR &_evaluator );
        //! Connects band-gap with quantities and printing policies
        template< class T_EVALUATOR > void bandgap( T_EVALUATOR &_evaluator );
        //! Connects band-gap with quantities and printing policies, and sets evaluator do_dipole.
        template< class T_EVALUATOR > void dipole( T_EVALUATOR &_evaluator );
        //! Connects electron mass with quantities and printing policies.
        template< class T_EVALUATOR > void emass( T_EVALUATOR &_evaluator );
        //! Connects hole mass with quantities and printing policies.
        template< class T_EVALUATOR > void hmass( T_EVALUATOR &_evaluator );
        //! Connects a property with quantities and printing policies.
        template< class T_EVALUATOR, class T_PROPERTY >
          details::ConnectProperty<T_EVALUATOR, T_PROPERTY > 
            connect_property( T_EVALUATOR &_evaluator, T_PROPERTY _property, 
                              const std::string &_name, const std::string &_decl );
        //! Reads objectives from XML input file.
        template< class T_FACTORY >
          void read_physical_properties( T_FACTORY& _factory, 
                                         boost::filesystem::path &_path );

        //! Prints out connections.
        void declare( const std::string& _string );


      } // namespace factory
    } // namespace Layered

  }
} // namespace LaDa

#include "factory.impl.hpp"

#endif
