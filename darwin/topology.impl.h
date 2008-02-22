//
//  Version: $Id$
//

namespace GA
{
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

  template< class T_CONTAINER >
    void Topology :: syncpop( T_CONTAINER &_cont )
    {
      return;
#ifdef _MPI
      //! In the case of a \e not \e graph, eg single pool topology, 
      //! all data should be equivalent across all processes at all times.
      //! the only exception is are printing and saving and such.
      if( not graph ) return;

      //! Farmhands never need do anything, except wait for the end of the run.
      if( graph->type == FARMHAND ) return;

      //! Cows receive individuals to analyze on a per-operation basis.
      //! They do not require full containers.
      if( graph->type == COW ) return;
      //! Bulls receive the container as broadcasted by the farmer.
      mpi::BroadCast bc( graph->farmer_comm );
      if( graph->type == BULL ) 
      {
        bc << mpi::BroadCast::allocate << mpi::BroadCast::broadcast;
        _cont.clear();
        typename T_CONTAINER :: value_type value;
        while( value->broadcast( bc ) ) _cont.push_back( value );
        return;
      }
      //! Finally, the farmer gets to broadcast the goods
      typename T_CONTAINER :: iterator i_var = _cont.begin();
      typename T_CONTAINER :: iterator i_var_end = _cont.end();
      for(; i_var != i_var_end; ++i_var )
        i_var->broadcast( bc );
      bc << mpi::BroadCast::allocate;
      for( i_var = _cont.begin(); i_var != i_var_end; ++i_var )
        i_var->broadcast( bc );
      bc << mpi::BroadCast::broadcast << mpi::BroadCast::reset();
#endif
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
  inline bool Topology :: objective () const
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
  inline bool Topology :: scaling() const
  {
    __SERIALCODE( return true; )
    __MPICODE( if( not graph ) return true;
               else if( graph->type == FARMER ) return true;
               return false;
             )
  }
  inline bool Topology :: store () const
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


  template <class T_GATRAITS> Breeder<T_GATRAITS>*
    Topology :: breeder ( typename T_GATRAITS :: t_Evaluator &_eval )
    {
      typedef T_GATRAITS t_GATraits;
      typedef typename t_GATraits :: t_Individual t_Individual;
                 
      __TRYCODE( __SERIALCODE( return new GA::Breeder<t_GATraits>; )
                 __MPICODE( if( not graph ) return new GA::Breeder<t_GATraits>();
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
    T_BASE<T_GATRAITS>* Topology :: evaluation ()
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
        t_History *result = NULL;
        __SERIALCODE( result = new t_History; )
        __TRYMPICODE( if( (not graph) ) result = new History< t_Individual >; 
                      else if( graph->type == FARMER ) 
                        result = new t_History;
                      else if( graph->type == BULL )
                        result = new mpi::Graph::BullHistory< t_GATraits >;,
                      "Error while creating history.\n" 
                    )
        if( result ) _eostates.storeFunctor( static_cast<t_History*>(result) );
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
    Topology :: special_store ( typename T_GATRAITS :: t_Evaluator& _eval )
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

