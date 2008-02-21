//
//  Version: $Id$
//
#ifndef  _SINGLE_TOPOLOGY_H_
#define  _SINGLE_TOPOLOGY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

#include <graph/topology.h>

namespace GA
{
  //! Wrapper around possible mpi topologies
  class Topology : private ::mpi::Base
  {
    protected:
    //! Seed of the random number generator;
    std::vector< types::t_int > seeds; 
   
    //! graph topology if created
    __MPICODE( Graph::Topology *graph; )

    public:  
      //! Constructor and Initializer.
      Topology() : ::mpi::Base::Base(), seeds(1,0)
                   __MPICODE(, graph(NULL) ) {}
      //! Copy Constructor.
      Topology   ( Topology &_comm ) : ::mpi::Base::Base(_comm) {}
      //! Destructor
      ~Topology();
  
      //! \brief Creates the GA mpi topology and groups.
      template < class T_EVALUATOR >
      bool Load( const TiXmlElement &_node, T_EVALUATOR &_eval );

      //! Loads seeds from attribute.
      bool LoadSeed( const TiXmlAttribute &_att );
      //! Print seeds to a string
      std::string print_seeds() const;
      //! Reseeds all procs as needed
      void reseed();
      //! Creates history objects
      template < class T_GATRAITS >
      History<typename T_GATRAITS::t_Individual> * history( eoStates & _eostates );
  };

  template< class T_EVALUATOR >
  bool Topology::Load( const TiXmlElement &_node, T_EVALUATOR &_eval )
  {
    __SERIALCODE( return true; )
    __MPICODE( graph = new Graph::Topology( _comm );
               if( graph->Load( _node ) and graph->init() ) return true;
               delete graph;
               graph = NULL;
             )
  }

  template< class T_GATRAITS > 
    History<typename T_GATRAITS :: t_Individual >*
      Topology :: history( eoStates & _eostates )
      {
        typedef T_GATRAITS t_GATraits;
        typedef t_GATraits :: t_Individual t_Individual;
        History<t_Individual> *result = NULL;
        __SERIALCODE( result = new History<t_Individual>; )
        __MPICODE( if( (not graph) ) result = new History< t_Individual >; 
                   else if( graph->type = FARMER ) result = new History< t_Individual >; 
                   else if( graph->type = BULL )
                     result = new mpi::Graph::BullHistory< t_GATraits >
                 )
        if( history ) _eostates.storeFunctor( history );
        return result;
      }
} // namespace GA


#endif
