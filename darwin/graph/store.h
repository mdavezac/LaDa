//
//  Version: $Id$
//
#ifndef  _DARWIN_COMMUNICATORS_H_
#define  _DARWIN_COMMUNICATORS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <darwin/store.h>

#ifdef _MPI

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
      class BullStore : protected Comm::Bull< T_GATRAITS, BullStore<T_GATRAITS> >,
                        public GA::Store::Base<T_GATRAITS>
      {
        public:
          typedef T_GATRAITS t_GATraits; //!< all GA classes \sa Traits::GA
        
        private:
          //! Type of the evaluator class .
          typedef typename t_GATraits :: t_Evaluator              t_Evaluator;
          //! Type of the base class.
          typedef GA::Store::Base<t_GATraits>                     t_Base;   
          //! Type of the Individual.
          typedef typename t_GATraits :: t_Individual             t_Individual;
          //! Type of the communication base class.
          typedef Comm::Bull< t_GATraits, BullStore<t_GATraits> > t_CommBase;

        public:
          //! Constructor and Initializer
          BullStore   (t_Evaluator &_eval, Graph::Topology *_topo )
                     :  t_CommBase( _topo ), t_Base( _eval) {}
          //! Destructor
          virtual ~BullStore() {}

          //! Transfers call to farmer.
          virtual void operator()( const t_Individual &_indiv );

          //! Does nothing
          bool Restart( const TiXmlElement &_node ) {}
          //! Does nothing
          bool Save( TiXmlElement &_node ) const {}

          //! Does nothing
          virtual void print_results(types::t_unsigned _age,
                                     bool is_comment = false) const {};
          //! Does nothing
          virtual std::string print() const { return ""; }
          //! Returns "GA::mpi::Graph::BullStore"
          virtual std::string what_is() const { return "GA::mpi::Graph::BullStore"; }
     
          //! Does nothing
          virtual void apply_all( eoMonOp<const t_Individual> *_op ) const {}
          //! Does nothing
          virtual void apply_best( eoMonOp<const t_Individual> *_op ) const {}
      };
      
      template < class T_GATRAITS >
        void BullStore<T_GATRAITS> :: operator()( const t_Individual &_indiv )
        {
          t_CommBase::request( t_CommBase::t_Requests::STORE );
          t_CommBase::send_individual( _indiv );
        }
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#endif
#endif
