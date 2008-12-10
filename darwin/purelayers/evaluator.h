//
//  Version: $Id$
//
#ifndef _DARWIN_PURELAYERS_EVALUATOR_H_
#define _DARWIN_PURELAYERS_EVALUATOR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>
#include <string>
#include <ostream>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include <crystal/structure.h>
#include <opt/types.h>
#include <mpi/mpi_object.h>
#include <vff/layered.h>

#include "../alloylayers/evaluator.h"

#include "translate.h"

namespace LaDa
{
  namespace GA
  {
    //! Optimization for superlattices of pure layers of quaternary alloys.
    namespace PureLayers
    {
      //! An evaluator class for pescan/vff in epitaxial alloy configuration space.
      template< class T_INDIVIDUAL,
                template<class> class T_TRANSLATE,
                template<class,class> class T_ASSIGN >
        class Evaluator : public AlloyLayers::Evaluator< T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN >
        {
            //! Type of the base class.
            typedef AlloyLayers::Evaluator< T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN > t_Base;
          public:
            //! Type of the individual.
            typedef typename t_Base :: t_Individual t_Individual;
            //! Type of the translation policy.
            typedef typename t_Base :: t_Translate t_Translate;
            //! Type of the assignement policy.
            typedef typename t_Base :: t_Assign t_Assign;
            //! All pertinent %GA traits
            typedef Traits::GA< Evaluator > t_GATraits;

            //! Concentration functor type.
            typedef GA::PureLayers::Concentration t_Concentration;
            //! Concentration functor.
            t_Concentration concentration;
       
            //! Constructor
            Evaluator() : t_Base() {}
            //! Copy Constructor
            Evaluator   ( const Evaluator &_c )
                      : t_Base(_c) {}
            //! Destructor
            virtual ~Evaluator() {};
       
            //! Loads structure, lattice, bandgap, vff from XML
            bool Load( const TiXmlElement &_node );
       
       
            using t_Base::Load;
            using t_Base::init;
            using t_Base::evaluate;
        };

    } // namespace Layered


 }
} // namespace LaDa
 

#include "evaluator.impl.h"

#endif 
