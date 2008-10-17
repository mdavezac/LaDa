//
//  Version: $Id$
//
#ifndef  _DARWIN_TOPOLOGY_H_
#define  _DARWIN_TOPOLOGY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/utils/eoState.h>
#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

#include "graph/graph.h"
#include "graph/objective.h"
#include "graph/evaluation.h"
#include "graph/history.h"
#include "graph/taboo.h"
#include "graph/store.h"
#include "graph/breeders.h"

#ifdef _MPI
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#endif

namespace GA
{
  //! Wrapper around possible mpi topologies
  class Topology
  {
    protected:
    //! Seed of the random number generator;
    std::vector< types::t_int > seeds; 
   
    //! graph topology if created
    __MPICODE( mpi::Graph::Topology *graph; )
    //! Head Communicator.
    __MPICODE( ::boost::mpi::communicator comm; )

    public:  
      //! Constructor and Initializer.
      Topology() : seeds(1,0)
                   __MPICONSTRUCTORCODE( graph(NULL) )
                   __MPICONSTRUCTORCODE( comm( *::mpi::main, 
                                               boost::mpi::comm_duplicate ) ) {}
      //! Copy Constructor.
      Topology   ( const Topology &_comm )
               : seeds( _comm.seeds )
                 __MPICONSTRUCTORCODE( graph(_comm.graph) ) 
                 __MPICONSTRUCTORCODE( comm( _comm.comm,
                                             boost::mpi::comm_attach) ) {}
#ifdef _MPI
      //! Constructor and Initializer.
      Topology   ( boost::mpi::communicator &_comm )
               : seeds(1,0), graph(NULL), 
                 comm( _comm, boost::mpi::comm_duplicate )  {}
#endif
      //! Destructor
      ~Topology() {};
  
      //! \brief Creates the GA mpi topology and groups.
      template < class T_EVALUATOR >
      bool Load( const TiXmlElement &_node, T_EVALUATOR &_eval );
      //! Print parameters of this run.
      std::string print() const;

      //! Loads seeds from attribute.
      bool LoadSeeds( const TiXmlAttribute &_att );
      //! Print seeds to a string
      std::string print_seeds() const;
      //! Reseeds all procs as needed
      void reseed();

      //! \brief Synchronizes populations across process which require it.
      //! \details Should also work for most containers, provided
      //!          T_CONTAINER::value_type::broadcast( mpi::BroadCast & )
      //!          exists.
      template< class T_CONTAINER > void syncpop( T_CONTAINER &_cont );

      //! Whether to create continuators (other than continuator container.)
      bool continuators() const;
      //! Whether to create a history.
      bool history() const;
      //! Whether to create mating operators.
      bool mating() const;
      //! Whether to create objectives.
      bool objective() const;
      //! Whether to populate.
      bool populate() const;
      //! Whether to create replacements.
      bool replacement() const;
      //! Whether to perform restart.
      bool restart() const;
      //! Whether to save from this proc.
      bool save() const;
      //! Whether to create scaling.
      bool scaling() const;
      //! Whether to create storage.
      bool store () const;
      //! Whether to create taboos
      bool taboos() const;

      //! returns breeder.
      template < class T_GATRAITS > GA::Breeder<T_GATRAITS>*
        breeder ( typename T_GATRAITS :: t_Evaluator &_eval );
      //! Creates an evaluation object.
      template <class T_GATRAITS, template <class> class T_BASE >  
        Evaluation :: Abstract< typename T_GATRAITS :: t_Population >*
          evaluation( typename T_GATRAITS :: t_Evaluator &_eval ); 
      //! Creates history objects
      template < class T_GATRAITS > History<typename T_GATRAITS::t_Individual>*
        history( eoState & _eostates );
      //! returns objective.
      template <class T_GATRAITS> typename GA::Objective::Types<T_GATRAITS>::t_Vector*
        objective ( const TiXmlElement &_node );

      //! Passes taboo pointer to farmer breeder.
      template <class T_GATRAITS> 
        void set ( GA::Breeder<T_GATRAITS> *_breeder,
                   Taboo_Base<typename T_GATRAITS::t_Individual> *_taboos);
      //! Creates special storage interfaces for very special people.
      template <class T_GATRAITS> typename GA::Store::Base<T_GATRAITS>*
        special_store ( typename T_GATRAITS :: t_Evaluator& _eval );
      //! Creates special_taboo interfaces for very special people.
      template < class T_GATRAITS > Taboo_Base<typename T_GATRAITS::t_Individual>*
        special_taboos( eoState & _e );
  };

} // namespace GA

#include "topology.impl.h"

#endif
