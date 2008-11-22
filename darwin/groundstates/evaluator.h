//
//  Version: $Id$
//
#ifndef _DARWIN_GROUNDSTATES_EVALUATOR_H_
#define _DARWIN_GROUNDSTATES_EVALUATOR_H_

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

#include "evaluator_base.h"
#include "../bitstring.h"
#include "../ce.h"

#include "translate.h"
#include "assign.h"

namespace LaDa
{
  namespace GA
  {
    namespace GroundStates
    {
      //! An evaluator class for cluster-expansion.
      template< class T_INDIVIDUAL,
                template<class> class T_TRANSLATE,
                template<class,class> class T_ASSIGN >
        class Evaluator : public EvaluatorBase< T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN >
        {
            //! Type of the base class.
            typedef EvaluatorBase< T_INDIVIDUAL, T_TRANSLATE, T_ASSIGN > t_Base;
          public:
            //! Type of the individual.
            typedef typename t_Base :: t_Individual t_Individual;
            //! Type of the translation policy.
            typedef typename t_Base :: t_Translate t_Translate;
            //! Type of the assignement policy.
            typedef typename t_Base :: t_Assign t_Assign;
            //! All pertinent %GA traits
            typedef Traits::GA< Evaluator > t_GATraits;
       
            //! Constructor
            Evaluator() : t_Base(), bandgap(structure) {}
            //! Copy Constructor
            Evaluator   ( const Evaluator &_c )
                      : t_Base(_c), bandgap(_c.bandgap) {}
            //! Destructor
            virtual ~Evaluator() {};
       
            //! Loads structure, lattice, bandgap, vff from XML
            bool Load( const TiXmlElement &_node );
       
            //! Initializes the bandgap.
            void init( t_Individual& _indiv ) { t_Base::init( _indiv ); ce.init(); }
       
            //! Calls ce evaluator.
            void evaluate();
             { t_Assign::assign( ce.evaluate(), current_individual->quantities() ); }

            //! Computes the gradient of the the functional
            void evaluate_gradient( t_QuantityGradients& _grad )
              { ce.evaluate_gradient( _grad ); }
            //! Evaluates the functional and computes the gradient
            void evaluate_with_gradient( t_QuantityGradients& _grad )
              { t_Assign::assign( ce.evaluate_with_gradient( _grad ), 
                                  current_individual->quantities() ); }
            //! \brief Computes component \a _pos of the gradient and stores in location \a
            //!        _pos of \a _grad
            void evaluate_one_gradient( t_QuantityGradients& _grad, types::t_unsigned _pos ) 
              { _grad[_pos] = ce.evaluate_one_gradient( _pos ); }
       
  #         ifdef _MPI
              //! forwards comm and suffix to bandgap.
              void set_mpi( boost::mpi::communicator *_comm, const std::string &_str )
                { ce.set_mpi( _comm, _str ); t_Base::set_mpi( _comm, _str ); } 
  #         endif
            using t_Base::Load;

            //! Loads a phenotypic niche from input
            //! \see  SingleSite::new_Niche_from_xml()
            void* Load_Niche( const TiXmlElement &_node )
              { return (void *) SingleSite::new_Niche_from_xml<t_GATraits, 1>
                                                              ( _node, concentration ); }

            //! Returns a constant reference to a structure.
            const Crystal::Structure& get_structure() const { return structure; }
       
          protected:
            //! Type of the cluster expansion interface.
            typedef CE::Darwin t_CE;
       
            //! Cluster expansion interface.
            t_CE ce; 
            
            using t_Base :: structure;
        };



    } // namespace Layered
    /** @} */


 }
} // namespace LaDa
 

#include "evaluator.impl.h"

#endif // _LAYERED_H_
