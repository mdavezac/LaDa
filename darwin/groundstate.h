//
//  Version: $Id$
//
#ifndef _GROUNDSTATE_H_
#define _GROUNDSTATE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <algorithm>
#include <functional>

#include <tinyxml/tinyxml.h>
#include <eo/eoOp.h>

#include <opt/types.h>
#include <mpi/mpi_object.h>

#include "evaluator.h"
#include "individual.h"
#include "single_site.h"
#include "ce.h"


/** \ingroup Genetic
 * @{ */
//! \brief Enables a single-objective search of a single-cell-shape
//!        cluster-expanded functional 
//! \details At this point, and as always, this one cell-shape  at a time. A
//!          physical individual is defined using the object, concentration, and
//!          fourier transforms of SingleSite. The evaluator class is mostly a
//!          wrapper around a CE::Darwin instance.
namespace GroundState
{
  using SingleSite::Object; 

  //! The physical individual type
  typedef Individual::Types< Object, 
                             SingleSite::Concentration, 
                             SingleSite::Fourier        > :: Scalar t_Individual;

  //! \brief Evaluator class for doing %GA over a cluster-expansion functional
  //! \xmlrestart the CE results are saved into an individual as the attribute CE.
  class Evaluator : public SingleSite::Evaluator< t_Individual >
  {
    public:
      //! Type of the individual
      typedef GroundState::t_Individual t_Individual;
      //! All %GA related types
      typedef Traits::GA< Evaluator > t_GATraits;
    protected:
      //! \cond
      typedef Evaluator t_This;
      typedef SingleSite::Evaluator< t_Individual > t_Base;
      //! \endcond
      
    public:
      using t_Base::Load;
      using t_Base::Save;

    protected:
      using t_Base :: current_individual;
      using t_Base :: current_object;

    protected:
      //! The interface to the single-cell-shape cluster expansion functional
      CE::Darwin ce; 

    public:
      //! Constructor
      Evaluator() : ce( structure ) {}; 
      //! Destructor
      ~Evaluator() {}

      //! Saves an individual to XML, including CE result
      bool Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      //! Loads an individual from XML, including CE result
      bool Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      //! Loads the functional from XML 
      bool Load( const TiXmlElement &_node );
      //! Loads a phenotypic niche from input
      //! \see  SingleSite::new_Niche_from_xml()
      void* Load_Niche( const TiXmlElement &_node )
        { return (void *) SingleSite::new_Niche_from_xml<t_GATraits, 1>( _node, concentration ); }

      //! Intializes before calls to evaluation member routines
      void init( t_Individual &_indiv )
        { t_Base :: init( _indiv ); ce.init( _indiv ); }

      //! \brief Evaluates the functional and stores result in the current
      //!        individual's quantity
      void evaluate()
        { current_individual->quantities() = ce.evaluate(); }
      //! Computes the gradient of the the functional
      void evaluate_gradient( t_QuantityGradients& _grad )
        { ce.evaluate_gradient( _grad ); }
      //! Evaluates the functional and computes the gradient
      void evaluate_with_gradient( t_QuantityGradients& _grad )
        { current_individual->quantities() = ce.evaluate_with_gradient( _grad ); }
      //! \brief Computes component \a _pos of the gradient and stores in location \a
      //!        _pos of \a _grad
      void evaluate_one_gradient( t_QuantityGradients& _grad, types::t_unsigned _pos ) 
        { _grad[_pos] = ce.evaluate_one_gradient( _pos ); }
  };

} // namespace GROUNDSTATE
/** @} */
#endif // _GROUNDSTATE_H_
