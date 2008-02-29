//
//  Version: $Id$
//
#ifndef  _GRAPH_HISTORY_H_
#define  _GRAPH_HISTORY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include <darwin/taboos.h>
#include "comm.h"

namespace GA
{
  namespace mpi 
  {
    namespace Graph
    {
      //! Transfer storage calls to the Farmer
      template< class T_GATRAITS >
      class BullHistory : protected Comm::Bull< T_GATRAITS, BullHistory<T_GATRAITS> >,
                          public ::GA::History<typename T_GATRAITS :: t_Individual>
      {
        public:
          typedef T_GATRAITS t_GATraits; //!< all GA classes \sa Traits::GA
        
        private:
          //! Type of the individual
          typedef typename t_GATraits :: t_Individual t_Individual;
          //! Type of the public base class.
          typedef ::GA::History< t_Individual >       t_Base;
          //! Type of this class.
          typedef BullHistory<t_GATraits>             t_This;
          //! Type of the communication base class.
          typedef Comm::Bull<t_GATraits, t_This>      t_CommBase;

        public:
          //! Constructor
          BullHistory( Topology *_topo ) : t_CommBase( _topo ), t_Base() {}
          //! Copy constructor
          BullHistory   ( const t_This &_taboo )
                    : t_Base(_taboo), t_CommBase( _taboo ){}
          //! Destructor
          virtual ~BullHistory() {};
    
          //! returns true if _indiv is in taboo_list
          virtual bool clone(t_Individual &_indiv );

          //! Should not be called as a taboo. throws.
          virtual bool operator()( const t_Individual& _indiv ) const
            { __THROW_ERROR( "Should not be called as a taboo.\n" ) }
      };
      
      template < class T_GATRAITS >
        bool BullHistory<T_GATRAITS> :: clone( t_Individual &_indiv )
        {
          t_CommBase::request( t_CommBase::t_Requests::TABOOCHECK );
          t_CommBase::send_individual( _indiv );
          bool result;
          t_CommBase::receive_object( result );
          if( not result ) return false;
          t_CommBase::receive_individual( _indiv );
          return result;
        }
    } // namespace Graph
  } // namespace mpi
} // namespace GA

#endif
#endif
