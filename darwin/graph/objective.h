//
//  Version: $Id$
//
#ifndef  _GRAPH_OBJECTIVE_H_
#define  _GRAPH_OBJECTIVE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <darwin/objective.h>

#ifdef _MPI

#include <list>
#include <utility>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

#include "comm.h"

namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      template< class T_GATRAITS >
      class BullObjective : protected Comm::Bull< T_GATRAITS,
                                                  BullObjective<T_GATRAITS> >,
                            public GA::Objective::Types<T_GATRAITS> :: t_Vector
      {
        public:
          //! all %GA traits
          typedef T_GATRAITS t_GATraits;
    
        protected:
          //! Type of individual in this %GA
          typedef typename t_GATraits :: t_Individual         t_Individual;
          //! Type of the fitness, as declared in the base class
          typedef typename t_GATraits :: t_Fitness            t_Fitness;
          //! Type of the quantity traits, as declared in the base class
          typedef typename t_GATraits :: t_QuantityTraits     t_QuantityTraits;
          //! Type of the quantity, as declared in the base class
          typedef typename t_QuantityTraits :: t_Quantity     t_Quantity;
          //! Type of the scalar quantity, as declared in the base class
          typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
          //! Type of the lamarckian traits, as declared in the base class
          typedef typename t_GATraits :: t_VA_Traits          t_VA_Traits;
          //! Type of the lamarckian gradients, as declared in the base class
          typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
          //! Type of the lamarckian variables, as declared in the base class
          typedef typename t_VA_Traits :: t_Type              t_VA_Type;
          //! A functor to for saving individuals
          typedef GA::SaveObject<t_GATraits>                  t_SaveOp;
          //! A functor to for loading individuals
          typedef GA::LoadObject<t_GATraits>                  t_LoadOp;
          //! Type of the communication base class.
          typedef Comm::Bull< t_GATraits, BullObjective<t_GATraits> > t_CommBase;
          //! Type of the objective base class.
          typedef typename GA::Objective::Types<t_GATraits> :: t_Vector t_Base;
        
        protected:
          using t_CommBase::TAG;
    
        public:
          //! Constructor.
          BullObjective   ( Topology *_topo )
                        : t_CommBase( _topo ) {}
          //! Copy Constructor.
          BullObjective   ( const BullObjective &_c )
                        : t_CommBase( _c ), t_Base( _c ) {}
    
          //! Calls upon each objective to save current status
          virtual bool Save( TiXmlElement &_node, t_SaveOp& _op) { return true; }
          //! Calls upon each objective to restart from a previously saved status
          virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op) { return true; }
          //! Returns true if at least one objective is impermanent
          virtual bool does_store() const { return false; }
          //! Returns "GA::mpi::Graph::Objective".
          virtual std::string print() const { return "GA::mpi::Graph::Objective"; }


          //! \brief Returns the fitness obtained from the farmer.
          virtual const t_Fitness& operator()(const t_Quantity& _q);
          //! \brief Returns the fitness obtained from the farmer.
          virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &,
                                                           t_QuantityGradients&,
                                                           t_VA_Type *);
          //! \brief Returns the fitness obtained from the farmer.
          virtual void evaluate_gradient( const t_Quantity &_q,
                                          t_QuantityGradients &_grad,
                                          t_VA_Type *_i_grad);
          //! \brief Returns the fitness obtained from the farmer.
          virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                                   t_QuantityGradients& _grad,
                                                   types::t_unsigned _n);
          //! Not needed by bull.
          virtual bool is_valid() const { return true; }
          //! Reeturns "GA::mpi::Graph::BullObjective".
          virtual std::string what_is() const 
            { return "GA::mpi::Graph::BullObjective"; }
      };
      
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#include "objective.impl.h"

#endif
#endif
