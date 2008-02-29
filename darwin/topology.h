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

#include "graph/topology.h"
#include "graph/objective.h"
#include "graph/evaluation.h"
#include "graph/history.h"
#include "graph/taboo.h"
#include "graph/store.h"
#include "graph/breeders.h"

namespace GA
{
  //! Wrapper around possible mpi topologies
  class Topology __MPICODE( : private ::mpi::Base )
  {
    protected:
    //! Seed of the random number generator;
    std::vector< types::t_int > seeds; 
   
    //! graph topology if created
    __MPICODE( mpi::Graph::Topology *graph; )

    public:  
#ifndef _MPI
      //! Constructor and Initializer.
      Topology() : seeds(1,0) {}
      //! Copy Constructor.
      Topology   ( const Topology &_comm ) : seeds( _comm.seeds ) {};
#else
      //! Constructor and Initializer.
      Topology() : ::mpi::Base::Base(), seeds(1,0),
                   graph(NULL)  {}
      //! Constructor and Initializer.
      Topology   ( ::mpi::Base &_comm )
               : ::mpi::Base::Base( _comm ), seeds(1,0),
                 graph(NULL)  {}
      //! Copy Constructor.
      Topology   ( const Topology &_comm ) 
               : ::mpi::Base::Base(_comm), seeds( _comm.seeds),
                 graph( _comm.graph ) {}
#endif
      //! Destructor
      ~Topology() {};
  
      //! \brief Creates the GA mpi topology and groups.
      template < class T_EVALUATOR >
      bool Load( const TiXmlElement &_node, T_EVALUATOR &_eval );

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
        T_BASE<T_GATRAITS>* evaluation (); 
      //! Creates history objects
      template < class T_GATRAITS > History<typename T_GATRAITS::t_Individual>*
        history( eoState & _eostates );
      //! returns objective.
      template <class T_GATRAITS> typename GA::Objective::Types<T_GATRAITS>::Vector*
        objective ( const TiXmlElement &_node );

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
