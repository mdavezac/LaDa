//
//  Version: $Id$
//
#ifndef  _SINGLE_TOPOLOGY_H_
#define  _SINGLE_TOPOLOGY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/utils/eoState.h>
#include <eo/eoBreed.h>
#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>

#include "graph/topology.h"
#include "objective.h"
#include "taboos.h"
#include "store.h"

namespace GA
{
  //! Wrapper around possible mpi topologies
  class Topology __MPICODE( : private ::mpi::Base )
  {
    protected:
    //! Seed of the random number generator;
    std::vector< types::t_int > seeds; 
   
    //! graph topology if created
    __MPICODE( Graph::Topology *graph; )

    public:  
#ifndef _MPI
      //! Constructor and Initializer.
      Topology() : seeds(1,0) {}
      //! Copy Constructor.
      Topology   ( Topology &_comm ) : seeds( _comm.seeds ) {};
#else
      //! Constructor and Initializer.
      Topology() : ::mpi::Base::Base(), seeds(1,0),
                   graph(NULL)  {}
      //! Copy Constructor.
      Topology   ( Topology &_comm ) 
               : ::mpi::Base::Base(_comm), seeds( _comm.seeds),
                 graph( _comm.graph ) {}
#endif
      //! Destructor
      ~Topology();
  
      //! \brief Creates the GA mpi topology and groups.
      template < class T_EVALUATOR >
      bool Load( const TiXmlElement &_node, T_EVALUATOR &_eval );

      //! Loads seeds from attribute.
      bool LoadSeeds( const TiXmlAttribute &_att );
      //! Print seeds to a string
      std::string print_seeds() const;
      //! Reseeds all procs as needed
      void reseed();

      __MPICODE( 
        //! \brief Synchronizes populations across process which require it.
        //! \details Should also work for most containers.
        template< T_CONTAINER > void syncpop( T_CONTAINER &_cont );
      )

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
      template < class T_GATRAITS > eoBreed<typename T_GATRAITS::t_Individual>*
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
        special_store ( typename T_GATRAITS :: t_Evaluator );
      //! Creates special_taboo interfaces for very special people.
      template < class T_GATRAITS > Taboo_Base<typename T_GATRAITS::t_Individual>*
        special_taboos( eoState & _e );
  };

  template< class T_EVALUATOR >
  bool Topology::Load( const TiXmlElement &_node, T_EVALUATOR &_eval )
  {
    __SERIALCODE( return true; )
    __TRYMPICODE( graph = new Graph::Topology( _comm );
                  if( graph->Load( _node ) and graph->init() ) return true;
                  delete graph;
                  graph = NULL;,
                  "Error while loading topology.\n"
                )
  }

  inline bool Topology :: history() const
  {
    __MPICODE( if(      graph and graph->type == COW ) return false;
               else if( graph and graph->type == FARMHAND ) return false;
             )
    return true;
  }
  inline bool Topology :: mating() const
  {
    __MPICODE( if(      graph and graph->type == FARMER ) return false;
               else if( graph and graph->type == COW ) return false;
               else if( graph and graph->type == FARMHAND ) return false;
             )
    return true;
  }
  bool Topology :: objective () const
  {
    __MPICODE( if( graph and graph->type == FARMHAND ) return false;
               if( graph and graph->type == COW ) return false;
             )
    return true;
  }
  inline bool Topology :: populate() const
  {
    __SERIALCODE( return true; )
    __MPICODE( if( not graph ) return true;
               else if( graph->type == FARMER ) return true;
               return false;
             )
  }
  inline bool Topology :: replacement() const
  {
    __SERIALCODE( return true; )
    __MPICODE( if( not graph ) return true;
               else if( graph->type == FARMER ) return true;
               return false;
             )
  }
  inline bool Topology :: restart() const
  {
    __SERIALCODE( return true; )
    __MPICODE( if( (not graph) and comm->rank() == 0 ) return true;
               else if( graph->type == FARMER ) return true;
               return false;
             )
  }
  inline bool Topology :: save() const
  {
    __SERIALCODE( return true; )
    __MPICODE( if( (not graph) and comm->rank() == 0 ) return true;
               else if( graph->type == FARMER ) return true;
               return false;
             )
  }
  bool Topology :: store () const
  {
    __MPICODE( if( graph and graph->type == FARMHAND ) return false;
               if( graph and graph->type == COW ) return false;
             )
    return true;
  }
  inline bool Topology :: taboos() const
  {
    __MPICODE( if( (not graph) ) return true;
               else if( graph->type == COW ) return false;
               else if( graph->type == FARMHAND` ) return false;
             )
    return true;
  }


  template <class T_GATRAITS> eoBreed<typename T_GATRAITS :: t_Individual> *
    breeder ( typename T_GATRAITS :: t_Evaluator &_eval )
    {
      typedef T_GATRAITS t_GATraits;
      typedef typename t_GATraits :: t_Individual t_Individual;
                 
      __TRYCODE( __SERIALCODE( return new GA::Breeder<t_Individual>(); )
                 __MPICODE( if( not graph ) return new GA::Breeder<t_Individual>();
                            if( graph->type == FARMER )
                              return new Breeders::Farm<t_GATraits>(graph);
                            if( graph->type == BULL )
                              return new Breeders::BULL<t_GATraits>(graph)
                            Breeders::Cow<t_GATraits > *breeder 
                              = Breeder::Cow<t_GATraits >(graph);
                            breeder->set( _eval );
                            return breeder;
                          ), 
                 "Error while creating Breeders.\n" 
               )
    }
  template <class T_GATRAITS, template <class> class T_BASE > 
    T_BASE<T_GATRAITS>* evaluation ()
    {
      typedef T_GATRAITS t_GATraits;
      typedef T_BASE<t_GATraits> t_Base;
                 
      try
      {
#ifndef _MPI
        return new t_Base();
#else
         if( not graph ) return new t_Base();
         if( graph->type == FARMER )
           return new Evaluation::Farm<t_GATraits, T_BASE>(graph);
         if( graph->type == BULL )
           return new Evaluation::Bull<t_GATraits, T_BASE>(graph);
         if( graph->type == FARMHAND )
           return new Evaluation::FarmHand<t_GATraits, T_BASE>(graph);
         Evaluation::Cow<t_GATraits, T_BASE> *evaluation; 
         evaluation = Evaluation::Cow<t_GATraits, T_BASE>(graph);
         evaluation->set( _eval );
         return evaluation;
#endif
      }
      __CATCHCODE(, "Error while creating Evaluation.\n" )
   }
  template< class T_GATRAITS > 
    History<typename T_GATRAITS :: t_Individual >*
      Topology :: history( eoState & _eostates )
      {
        typedef T_GATRAITS t_GATraits;
        typedef typename t_GATraits :: t_Individual t_Individual;
        typedef History<t_Individual> t_History;
        History<t_Individual> *result = NULL;
        __SERIALCODE( result = new t_History; )
        __TRYMPICODE( if( (not graph) ) result = new History< t_Individual >; 
                      else if( graph->type == FARMER ) 
                        result = new t_History;
                      else if( graph->type == BULL )
                        result = new mpi::Graph::BullHistory< t_GATraits >;,
                      "Error while creating history.\n" 
                    )
        if( history ) _eostates.storeFunctor( static_cast<t_History*>(history) );
        return result;
      }
  template <class T_GATRAITS> typename GA::Objective::Types<T_GATRAITS>::Vector*
    Topology :: objective ( const TiXmlElement &_node )
    {
      typedef typename GA::Objective::Types<T_GATRAITS> t_ObjectiveType;
      __TRYCODE( 
        __MPICODE( if (      graph and graph->type == COW ) return NULL;
                   else if ( graph and graph->type == FARMHAND ) return NULL;
                   else if ( graph and graph->type == BULL )
                       return new BullObjective( graph ); 
                 )
        const TiXmlElement *child = _node.FirstChildElement("Objective");
        if ( not child ) child = _node.FirstChildElement("Method");
        return t_ObjectiveType :: new_from_xml( *child );,
        " Could not find Objective tag in input file.\n" 
      )
    }

  template <class T_GATRAITS> typename GA::Store::Base<T_GATRAITS>*
    Topology :: special_store ( typename T_GATRAITS :: t_Evaluator )
    {
      __TRYMPICODE( if( graph and graph->type == BULL ) 
                      return new BullStore( graph, _eval );,
                    "Error while creating BullStore.\n"
                  )
      return NULL;
    }
  template< class T_GATRAITS > 
    Taboo_Base<typename T_GATRAITS :: t_Individual >*
      Topology :: special_taboos( eoState &_e )
      {
        typedef T_GATRAITS t_GATraits;
        typedef typename t_GATraits :: t_Individual t_Individual;
        __SERIALCODE( return NULL; )
        __TRYMPICODE( if( (not graph) ) return NULL;
                      else if( graph->type == FARMER ) return NULL;
                      else if( graph->type == COW ) return NULL;
                      else if( graph->type == FARMHAND ) return NULL;
                      Taboo_Base<t_Individual> taboo 
                        = new mpi::Graph::BullTaboo< t_GATraits >;
                      _e.storeFunctor(taboo);,
                      "Error while creating taboos.\n"
                    )
      }

} // namespace GA


#endif
