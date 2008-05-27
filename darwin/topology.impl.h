//
//  Version: $Id$
//

namespace GA
{
  template< class T_EVALUATOR >
  bool Topology::Load( const TiXmlElement &_node, T_EVALUATOR &_eval )
  {
    __SERIALCODE( return true; )
    __TRYMPICODE(
      graph = new mpi::Graph::Topology( comm );
      if( graph->Load( _node ) )
        if( graph->init() )
        { 
          graph->set_mpi( _eval );
          return true;
        }
      delete graph;
      graph = NULL;
      std::string str = "";
      _eval.set_mpi( &comm, str );,
      "Error while loading topology.\n"
    )
  }

  template< class T_CONTAINER >
    void Topology :: syncpop( T_CONTAINER &_cont )
    {
#ifdef _MPI
      //! In the case of a \e not \e graph, eg single pool topology, 
      //! all data should be equivalent across all processes at all times.
      //! the only exception are printing and saving and such.
      if( not graph ) return;

      //! Farmhands never need do anything, except wait for the end of the run.
      if( graph->type == mpi::Graph::t_Type::FARMHAND ) return;

      //! Cows receive individuals to analyze on a per-operation basis.
      //! They do not require full containers.
      if( graph->type == mpi::Graph::t_Type::COW ) return;
      //! Broadcast from farmer to bulls
      Print :: out << "Broadcasting population: " << _cont.size() << Print::endl;
      boost::mpi::broadcast( *graph->farmer_comm(), _cont, 0 );
      Print :: out << "Broadcasted population: " << _cont.size() << Print::endl;
#endif
    }

  inline bool Topology :: continuators() const
  {
    __MPICODE(
      if( graph and graph->type != mpi::Graph::t_Type::FARMER ) return false;
    )
    return true;
  }
  inline bool Topology :: history() const
  {
    __MPICODE(
      if(      graph and graph->type == mpi::Graph::t_Type::COW ) return false;
      else if( graph and graph->type == mpi::Graph::t_Type::FARMHAND ) return false;
    )
    return true;
  }
  inline bool Topology :: mating() const
  {
    __MPICODE(
      if(      graph and graph->type == mpi::Graph::t_Type::FARMER ) return false;
      else if( graph and graph->type == mpi::Graph::t_Type::COW ) return false;
      else if( graph and graph->type == mpi::Graph::t_Type::FARMHAND ) return false;
    )
    return true;
  }
  inline bool Topology :: objective () const
  {
    __MPICODE(
      if( graph and graph->type == mpi::Graph::t_Type::FARMHAND ) return false;
      if( graph and graph->type == mpi::Graph::t_Type::COW ) return false;
    )
    return true;
  }
  inline bool Topology :: populate() const
  {
    __SERIALCODE( return true; )
    __MPICODE(
      if( not graph ) return true;
      else if( graph->type == mpi::Graph::t_Type::FARMER ) return true;
      return false;
    )
  }
  inline bool Topology :: replacement() const
  {
    __SERIALCODE( return true; )
    __MPICODE(
      if( not graph ) return true;
      else if( graph->type == mpi::Graph::t_Type::FARMER ) return true;
      return false;
    )
  }
  inline bool Topology :: restart() const
  {
    __SERIALCODE( return true; )
    __MPICODE(
      if( (not graph) and comm.rank() == 0 ) return true;
      else if( graph->type == mpi::Graph::t_Type::FARMER ) return true;
      return false;
    )
  }
  inline bool Topology :: save() const
  {
    __SERIALCODE( return true; )
    __MPICODE(
      if( (not graph) and comm.rank() == 0 ) return true;
      else if( graph->type == mpi::Graph::t_Type::FARMER ) return true;
      return false;
    )
  }
  inline bool Topology :: scaling() const
  {
    __SERIALCODE( return true; )
    __MPICODE(
      if( not graph ) return true;
      else if( graph->type == mpi::Graph::t_Type::FARMER ) return true;
      return false;
    )
  }
  inline bool Topology :: store () const
  {
    __MPICODE(
      if( graph and graph->type == mpi::Graph::t_Type::FARMHAND ) return false;
      if( graph and graph->type == mpi::Graph::t_Type::COW ) return false;
    )
    return true;
  }
  inline bool Topology :: taboos() const
  {
    __MPICODE(
      if( (not graph) ) return true;
      else if( graph->type == mpi::Graph::t_Type::COW ) return false;
      else if( graph->type == mpi::Graph::t_Type::FARMHAND ) return false;
    )
    return true;
  }


  template <class T_GATRAITS> GA::Breeder<T_GATRAITS>*
    Topology :: breeder ( typename T_GATRAITS :: t_Evaluator &_eval )
    {
      typedef T_GATRAITS t_GATraits;
      typedef typename t_GATraits :: t_Individual t_Individual;
      typedef GA::Breeder<T_GATRAITS> t_Return;
                 
      __TRYCODE(
        __SERIALCODE( return new GA::Breeder<t_GATraits>; )
        __MPICODE(
           if( not graph ) return new GA::Breeder<t_GATraits>();
           if( graph->type == mpi::Graph::t_Type::FARMER )
             return (t_Return*) new mpi::Graph::Breeders::Farmer<t_GATraits>(graph);
           if( graph->type == mpi::Graph::t_Type::BULL )
             return (t_Return*) new mpi::Graph::Breeders::Bull<t_GATraits>(graph);
           if( graph->type == mpi::Graph::t_Type::FARMHAND )
             return (t_Return*) new mpi::Graph::Breeders::Farmhand<t_GATraits>;
           mpi::Graph::Breeders::Cow<t_GATraits> *breeder 
               = new mpi::Graph::Breeders::Cow<t_GATraits>(graph);
           breeder->set( _eval );
           return (t_Return*) breeder;
         ), 
         "Error while creating Breeders.\n" 
       )
    }
  template <class T_GATRAITS, template <class> class T_BASE >  
    T_BASE<T_GATRAITS>* 
      Topology :: evaluation ( typename T_GATRAITS :: t_Evaluator &_eval )
      {
        typedef T_GATRAITS t_GATraits;
        typedef T_BASE<t_GATraits> t_Base;
#ifdef _MPI
        typedef mpi::Graph::Evaluation::Farmer<t_GATraits, T_BASE> t_Farmer;
        typedef mpi::Graph::Evaluation::Bull<t_GATraits, T_BASE> t_Bull;
        typedef mpi::Graph::Evaluation::Farmhand< t_Base > t_Farmhand;
        typedef mpi::Graph::Evaluation::Cow<t_GATraits, T_BASE> t_Cow;
#endif
//
//       __TRYCODE(
//         __SERIALCODE( return new t_Base(); )
//         __MPICODE( 
            if( not graph )
            {
              Print :: out << "topology.evaluation " << Print::endl;
              t_Base *result = new t_Base;
              result->set( _eval );
              return result;
            }
            if( graph->type == mpi::Graph::t_Type::FARMER )
            {
              Print :: out << "topology.evaluation Farmer "  << Print::endl;
              t_Farmer *result = new t_Farmer( graph );
              result->set( _eval );
              return (t_Base*) result;
            }
            if( graph->type == mpi::Graph::t_Type::BULL )
            {
              Print :: out << "topology.evaluation Bull "  << Print::endl;
              t_Bull *result = new t_Bull( graph );
              result->set( _eval );
              return (t_Base*) result;
            }
            if( graph->type == mpi::Graph::t_Type::FARMHAND )
            {
              Print :: out << "topology.evaluation Farmhand "  << Print::endl;
              t_Farmhand *result = new t_Farmhand;
              result->set( _eval );
              return (t_Base*) result;
            }
            t_Cow *result = new t_Cow( graph );
            result->set( _eval );
            return (t_Base*) result;
            Print :: out << "topology.evaluation Cow "  << Print::endl;
 //        ),
 //        "Error while creating Evaluation.\n" 
 //      )
     }
  template< class T_GATRAITS > 
    History<typename T_GATRAITS :: t_Individual >*
      Topology :: history( eoState & _eostates )
      {
        typedef T_GATRAITS t_GATraits;
        typedef typename t_GATraits :: t_Individual t_Individual;
        typedef History<t_Individual> t_History;
        t_History *result = NULL;
        __SERIALCODE( result = new t_History; )
        __TRYMPICODE(
          if( (not graph) ) result = new t_History; 
          else if( graph->type == mpi::Graph::t_Type::FARMER ) 
            result = new t_History;
          else if( graph->type == mpi::Graph::t_Type::BULL )
            result = (t_History*) new mpi::Graph::BullHistory< t_GATraits >(graph);,
          "Error while creating history.\n" 
        )
        if( result ) _eostates.storeFunctor( static_cast<t_History*>(result) );
        return result;
      }
  template <class T_GATRAITS> typename GA::Objective::Types<T_GATRAITS>::Vector*
    Topology :: objective ( const TiXmlElement &_node )
    {
      typedef T_GATRAITS t_GATraits;
      typedef typename GA::Objective::Types<t_GATraits> t_ObjectiveType;
      __TRYCODE( 
        __MPICODE(
          if (      graph and graph->type == mpi::Graph::t_Type::COW ) return NULL;
          else if ( graph and graph->type == mpi::Graph::t_Type::FARMHAND )
            return NULL;
          else if ( graph and graph->type == mpi::Graph::t_Type::BULL )
            return new mpi::Graph::BullObjective<t_GATraits>( graph ); 
        )
        const TiXmlElement *child = _node.FirstChildElement("Objective");
        if ( not child ) child = _node.FirstChildElement("Method");
        return t_ObjectiveType :: new_from_xml( *child );,
        " Could not find Objective tag in input file.\n" 
      )
    }

  template <class T_GATRAITS> typename GA::Store::Base<T_GATRAITS>*
    Topology :: special_store ( typename T_GATRAITS :: t_Evaluator& _eval )
    {
      __TRYMPICODE(
        if( graph and graph->type == mpi::Graph::t_Type::BULL ) 
          return new mpi::Graph::BullStore<T_GATRAITS>( _eval, graph );,
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
        __TRYMPICODE(
          if( (not graph) ) return NULL;
          else if( graph->type == mpi::Graph::t_Type::FARMER )
            return NULL;
          else if( graph->type == mpi::Graph::t_Type::COW )
            return NULL;
          else if( graph->type == mpi::Graph::t_Type::FARMHAND )
            return NULL;
          mpi::Graph::BullTaboo< t_GATraits > *taboo 
            = new mpi::Graph::BullTaboo< t_GATraits >(graph);
          _e.storeFunctor(taboo);
          return (Taboo_Base<t_Individual>*) taboo;,
          "Error while creating taboos.\n"
        )
      }

} // namespace GA

