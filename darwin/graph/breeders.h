//
//  Version: $Id$
//
#ifndef  _GRAPH_BRREEDERS_H_
#define  _GRAPH_BRREEDERS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <darwin/breeder.h>

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
      //! Contains all Breeder related stuff in the mpi::Graph Topology.
      namespace Breeders
      {

        //! A breeder class which does nothing.
        template<class T_GATRAITS>
        class Farmhand : public GA::Breeder<T_GATRAITS>
        {
          public:
            typedef T_GATRAITS t_GATraits; //!< all %GA traits
    
          protected:
            //! all individual traits
            typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
            //! type of an individual
            typedef typename t_GATraits::t_Individual  t_Individual; 
            //! type of the population
            typedef typename t_GATraits::t_Population  t_Population; 
            //! Base class type.
            typedef GA::Breeder<t_GATraits> t_Base;
    
          public:
            //! Constructor.
            Farmhand() : t_Base() {};
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents,
                            t_Population& _offspring) {};
       
            ///! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Farmhand"; }
        };
    
        //! \brief A breeder class to rule them all.
        //! \details This functor dispatches commands to the bulls, such as breed
        //!          one and stop breeding. 
        template<class T_GATRAITS>
        class Farmer : protected Comm::Farmer< T_GATRAITS, Farmer<T_GATRAITS> >, 
                       public GA::Breeder<T_GATRAITS>
        {
          friend class Comm::Farmer< T_GATRAITS, Farmer<T_GATRAITS> >;
          public:
            typedef T_GATRAITS t_GATraits; //!< all %GA traits
    
          protected:
            //! all individual traits
            typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
            //! type of an individual
            typedef typename t_GATraits::t_Individual  t_Individual; 
            //! type of the population
            typedef typename t_GATraits::t_Population  t_Population; 
            //! Base class type.
            typedef GA::Breeder<t_GATraits> t_Base;
            //! Communication base class
            typedef Comm::Farmer< T_GATRAITS, Farmer<t_GATraits> > t_CommBase;
    
          protected:
            types::t_unsigned target;
            t_Population *offspring;
            
          protected:
            using t_CommBase :: TAG;
    
          public:
            //! Constructor.
            Farmer   ( Topology *_topo )
                   : t_CommBase( _topo ), t_Base(),
                     target(0), offspring(NULL) {};
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Farmer"; }

            //! Sets taboo pointer
            void set( Taboo_Base<t_Individual> *_taboo )
              { t_CommBase::set( _taboo ); }
    
            //! \brief Response to WAITING request.
            //! \details Receives the individual. If still more are to be
            //!          created, sends t_CommBase::t_Commands::DONE, otherwise
            //!          t_CommBase::t_Commands::GO. The persistent request is
            //!          only activated in the latter case.
            void onWait( types::t_int _bull );
        };
    
        template<class T_GATRAITS>
        class Bull : protected Comm::Bull< T_GATRAITS, Bull<T_GATRAITS> >,
                     public GA::Breeder<T_GATRAITS>
        {
          friend class Comm::Bull< T_GATRAITS, Bull<T_GATRAITS> >;
          public:
            typedef T_GATRAITS t_GATraits; //!< all %GA traits
    
          protected:
            //! all individual traits
            typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
            //! type of an individual
            typedef typename t_GATraits::t_Individual  t_Individual; 
            //! type of the population
            typedef typename t_GATraits::t_Population  t_Population; 
            //! Base class type.
            typedef GA::Breeder<t_GATraits> t_Base;
            //! Communication base class
            typedef Comm::Bull< T_GATRAITS, Bull<t_GATraits> > t_CommBase;
    
          protected:
            using t_CommBase :: TAG;
    
          public:
            //! Constructor.
            Bull   ( Topology *_topo )
                 : t_CommBase(_topo ), t_Base() {}
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Breeder::Bull"; }
        };
    
        //! Just a typedef for Comm::BaseCow.
        template< class T_GATRAITS >
        class Cow : public Comm::LaNormande< T_GATRAITS, GA::Breeder > 
        {
          public:  
            //! All %GA types.
            typedef T_GATRAITS t_GATraits;
            //! Constructor.
            Cow   ( Topology * _topo )
                : Comm::LaNormande< t_GATraits, GA::Breeder >( _topo ){}
            using Comm::LaNormande< T_GATRAITS, GA::Breeder > :: set;
        };
      } // namespace Breeders
    } // namespace Graph
    
  } // namespace mpi
} // namespace GA

#include "breeders.impl.h"

#endif
#endif
