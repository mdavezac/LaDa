//
//  Version: $Id$
//
#ifndef  _GRAPH_TABOO_H_
#define  _GRAPH_TABOO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef _MPI

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include "comm.h"

namespace LaDa
{
  namespace GA
  {
    namespace mpi 
    {
      namespace Graph
      {
        //! Transfer storage calls to the Farmer
        template< class T_GATRAITS >
        class BullTaboo : protected Comm::Bull< T_GATRAITS, BullTaboo<T_GATRAITS> >
        {
          public:
            typedef T_GATRAITS t_GATraits; //!< all GA classes \sa Traits::GA
          
          private:
            //! Type of the individual
            typedef typename t_GATraits :: t_Individual t_Individual;
            //! Type of this class.
            typedef BullTaboo<t_GATraits> t_This;
            //! Type of the communication base class.
            typedef Comm::Bull< t_GATraits, t_This > t_CommBase;

          protected:
            using t_CommBase :: TAG; 

          public:
            //! Constructor
            BullTaboo ( Topology* _topo ) : t_CommBase(_topo) {}
            //! Copy constructor
            BullTaboo   ( const t_This &_taboo )
                      : t_CommBase( _taboo ){}
            //! Destructor
            virtual ~BullTaboo() {};
      
            //! returns true if _indiv is in taboo_list.
            bool operator()( const t_Individual& _indiv ) const; 
        };
        
        template < class T_GATRAITS >
          bool BullTaboo<T_GATRAITS> :: operator()( const t_Individual &_indiv ) const
          {
            t_CommBase::request( t_CommBase::t_Requests::TABOOCHECK );
            t_CommBase::comm->send( 0, ONTABOO_TAG1( TAG ), _indiv );
            bool result;
            t_CommBase::comm->recv( 0, ONTABOO_TAG2( TAG ), result );
            return result;
          }
      } // namespace Graph
    } // namespace mpi
  } // namespace GA
} // namespace LaDa
#endif
#endif
