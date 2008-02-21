//
//  Version: $Id$
//
#ifndef _GRAPH_COMM_H_
#define _GRAPH_COMM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef _MPI

#include <opt/types.h>
#include <opt/debug.h>
#include <mpi/mpi_object.h>
#include "topology.h"

namespace GA
{

  //! MPI for the Genetic Algorithm.
  namespace mpi
  {
    namespace Graph
    {
      //! \brief Holds all things communication related for the GA::mpi::Graph
      //!        topology
      //! \details In practice, it includes a number of helper functions which
      //!          can communicate object point to point. More importantly, it
      //!          includes base class meant to be used in CRT format
      //!          (Curriously Recurring Templates) for the breeder classes,
      //!          the objective classes, the taboo classes, the evaluation
      //!          classes, and the history classes. 
      namespace Comm
      {
        //! Sends a template object to bull \a _bull.
        template< class T_QUANTITY >
        void send_quantity( types::t_unsigned _bull,
                            T_QUANTITY &_q, MPI::Comm *_comm);
        //! Receives a template object from bull \a _bull.
        template< class T_QUANTITY >
        void receive_quantity( types::t_unsigned _bull,
                               T_QUANTITY &_q, MPI::Comm *_comm);
        //! Sends a template object to bull \a _bull.
        template< class T_OBJECT >
        void send_template_object( types::t_unsigned _bull,
                                   T_OBJECT &_o, MPI::Comm *_comm);
        //! Receives a template object from bull \a _bull.
        template< class T_OBJECT >
        void receive_template_object( types::t_unsigned _bull,
                                      T_OBJECT &_o, MPI::Comm *_comm);
        //! Sends a (non-template) object to bull \a _bull.
        template< class T_OBJECT >
        void send_object( types::t_unsigned _bull,
                          T_OBJECT _o, MPI::Comm *_comm);
        //! Receives a (non-template) object from bull \a _bull.
        template< class T_OBJECT >
        void receive_object( types::t_unsigned _bull,
                             T_OBJECT &_o, MPI::Comm *_comm);
        //! Sends a range (must be iterators or pointers).
        template< class T_ITERATOR >
        void send_range( types::t_unsigned _bull, T_ITERATOR _first,
                         T_ITERATOR  _last, MPI::Comm *_comm);
        //! Receive a range (must be iterators or pointers).
        template< class T_ITERATOR >
        void receive_range( types::t_unsigned _bull, T_ITERATOR _first,
                            T_ITERATOR  _last, MPI::Comm *_comm);
        //! Broadcasts a template object from \a _root.
        template< class T_OBJECT >  
        void bcast_template_object( types::t_int _root,
                                    T_OBJECT &_object, MPI::Comm *_comm )
        //! Broadcasts a range of iterators/pointers.
        template< class T_ITERATOR >  
          void bcast_range( types::t_int _root, T_ITERATOR _first,
                            T_ITERATOR _last, MPI::Comm *_comm )
        //! Broadcasts a template object from \a _root.
        template< class T_OBJECT >  
        void bcast_object( T_OBJECT &_object, MPI::Comm *_comm )
        
        //! \cond
        template< T_GATRAITS, T_DERIVED > class Cow;
        template< T_GATRAITS, T_DERIVED > class Bull;
        //! \endcond
        
        //! \brief Communication CRT class for farmers.
        //! \details This class defines a number of persistent requests a farmer
        //!          will expect to receive from its bulls. More specifically,
        //!          he will expect:
        //!          - Farmer::t_Requests::WAITING for when a
        //!            bull does not know what it will do next
        //!          - Farmer::t_Requests::OBJECTIVE for when a
        //!            bull needs the evaluation of an objective.
        //!          - Farmer::t_Requests::OBJECTIVE_GRADIENT for when a
        //!            bull needs the evaluation of the gradient of an
        //!            objective.
        //!          - Farmer::t_Requests::OBJECTIVE_WITH_GRADIENT for when a
        //!            bull needs the evaluation of an objective and its the
        //!            gradient.
        //!          - Farmer::t_Requests::OBJECTIVE_ONE_GRADIENT for when a
        //!            bull needs the evaluation of the gradient of an
        //!            objective in one particular direction.
        //!          - Farmer::t_Requests::TABOOCHECK for when a
        //!            bull needs to know whether a specific individual is
        //!            taboo or not.
        //!          - Farmer::t_Requests::HISTORYCHECK for when a
        //!            bull needs to know whether a specific individual is
        //!            already known or not.
        //!          - Farmer::t_Requests::STORE for a bull to propose an
        //!            individual for storage
        //!          .
        //!          The CRT will make calls to onWaiting(), onTaboo(),
        //!          onObjective(), onHistory() for each of these requests
        //!          respectively. These routines should be defined in the
        //!          derived classes. Defaults are provided for the last three.
        template< T_DERIVED >
        class Farmer : private ::mpi::Base
        {
          friend class Bull<T_GATRAITS, T_DERIVED>;
          friend class Cow<T_GATRAITS, T_DERIVED>;
          public:
            //! Type of the derived class
            typedef T_DERIVED t_Derived;
            //! All %GA traits
            typedef typename T_GATRAITS                           t_GATraits;
          private:
            //! Type of this class
            typedef Farmer<t_GATraits>                            t_This;
            //! Type of the individuals
            typedef typename t_GATraits :: t_Individual           t_Individual;
            //! %Traits of the quantity (or raw fitness)
            typedef typename t_GATraits :: t_QuantityTraits       t_QuantityTraits;
            //! Type of the genomic object
            typedef typename t_GATraits :: t_Object               t_Object;
            //! Type of the population
            typedef typename t_GATraits :: t_Population           t_Population;
            //! Type of the collection of populations
            typedef typename t_GATraits :: t_Islands              t_Islands;
            //! Type of the scalar quantity (or scalar raw fitness)
            typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
            //! Type of the objective type holder
            typedef typename Objective :: Types < t_GATraits >    t_ObjectiveType;
            //! Type of the storage type holder
            typedef typename Store :: Types< t_GATraits >         t_Store;
            //! Type of the fitness object.
            typedef typename t_IndivTraits :: t_Fitness           t_Fitness;
            //! Type of the base class
            typedef ::mpi::Base t_Base;
          protected:
            //! Possible requests from bulls.
            enum t_Requests
            {
              WAITING, //!< Waiting for what to do next.
              OBJECTIVE, //!< Requesting an objective evaluation.
              OBJECTIVE_GRADIENT, //!< Requesting an objective gradient evaluation.
              OBJECTIVE_WITH_GRADIENT, //!< Requesting an objective evaluation with gradient.
              OBJECTIVE_ONE_GRADIENT, //!< Requesting an objective evaluation of one gradient.
              TABOOCHECK, //!< Requesting to know whether an individual is taboo.
              HISTORYCHECK, //!< Requesting to know whether an individual is known.
              STORE; //!< Proposes an individual for storage.
            };
            enum t_Commands
            {
              GO, //!< Keep going with loop.
              DONE; //!< Exit loop.
            };
            const static types::t_int TAG = 1;
          
          protected:
            //! Holds requests from bull;
            MPI::Prequest *from_bulls;
            //! Number of bulls.
            types::t_unsigned nbulls;
            //! Request buffers
            types::t_unsigned *requests;
            //! Taboo functor.
            Taboo_Base<t_Individual>*          taboos;
            //! Objective functor.
            typename t_ObjectiveType::Vector*  objective;
            //! Store functor.
            typename t_Store :: Base*          store;
            //! History functor
            History<t_Individual>*             history;
        
          protected:
            Farmer( MPI :: Comm *_comm );
            ~Farmer();
        
            //! \brief Tests incoming messages from bull and dispatches response.
            //! \details The incoming messages are listed in
            //!          Farmer::t_Requests. There response are dispatched to
            //!          to onWait(), onTaboo(), onHistory(), and onCheck() through
            //!          a CRT mechanism (Curriuosly Recurring Template), e.g. to
            //!          the member functions of the derivatives of this class.
            bool test_bulls();
            //! Wait for bulls to finish.
            bool wait_bulls();
        
            //! Sends an individual to bull \a _bull.
            void send_individual( types::t_unsigned _bull, t_Individual &_indiv)
              { send_template_object< t_Individual >( _bull + 1, _indiv, 
                                                      t_CommBase::comm ); }
            //! Sends an individual to bull \a _bull.
            void receive_individual( types::t_unsigned _bull, t_Individual &_indiv )
              { receive_template_object< t_Individual >( _bull + 1, _indiv,
                                                         t_CommBase::comm ); }
            //! Sends an individual to bull \a _bull.
            void send_fitness( types::t_unsigned _bull, t_Fitness &_fit )
              { send_template_object< t_Fitness >( _bull + 1, _fit,
                                                   t_CommBase::comm ); }
            //! Sends an individual to bull \a _bull.
            void receive_individual( types::t_unsigned _bull, t_Fitness &_fit )
              { receive_template_object< t_Fitness >( _bull + 1, _fit,
                                                      t_CommBase::comm ); }
            //! Sends an individual to bull \a _bull.
            void send_quantities( types::t_unsigned _bull, t_Quantities &_q )
              { send_quantities< t_Fitness >( _bull + 1, _q, t_CommBase::comm ); }
            //! Sends an individual to bull \a _bull.
            void receive_quantities( types::t_unsigned _bull, t_Quantities &_fit )
              { receive_quantities< t_Fitness >( _bull + 1, _q, t_CommBase::comm ); }
            //! Sends gradients to bull.
            void send_gradients( types::t_unsigned _bull, t_VA_type *_grad,
                                 types::t_unsigned _n  );
              { send_range< t_VA_type >( _bull, _grad, _grad + _n, t_CommBase::comm ); }
            //! Receives gradients from bull.
            void receive_gradients( types::t_unsigned _bull, t_QuantityGradients &_grad) 
              { receive_quantities< t_QuantityGradients >( _bull, _grad, t_CommBase::comm ); }
            //! Sends an object to \a _bull.
            template < class t_OBJECT > void send_object( types::t_unsigned _bull,
                                                          T_OBJECT _object );
              { send_object< T_OBJECT >( _bull, _object, t_CommBase::comm ); }
            //! Receives an object from \a _bull.
            template < class t_OBJECT > void receive_object( types::t_object _bull, 
                                                             T_OBJECT &_object );
              { receive_object< T_OBJECT >( _bull, _object, t_CommBase::comm ); }
            //! Starts all persistent requests from bulls ( Farmer::requests )
            void startall() { _comm->StartAll( nbulls, requests ); } 
            //! Sends a command to \a _bull.
            void send_command( types::t_unsigned _bull, t_Commands _c );
            //! Activates request for \a _bull.
            void activate( types::t_unsigned _bull)
              { request[_bull].start(); }
            
            //! Response to TABOOCHECK request
            void onTaboo( types::t_int _bull );
            //! Response to OBJECTIVE request
            void onObjective( types::t_int _bull );
            //! Response to OBJECTIVE_GRADIENT request
            void onGradient( types::t_int _bull );
            //! Response to OBJECTIVE_WITH_GRADIENT request
            void onWithGradient( types::t_int _bull );
            //! Response to OBJECTIVE_ONE_GRADIENT request
            void onOneGradient( types::t_int _bull );
            //! Response to HISTORYCHECK request
            void onHistory( types::t_int _bull );
            //! Response to STORE request
            void onStore( types::t_int _bull );
        
          public:
            //! Sets taboo pointer
            set( Taboo_Base<t_Individual> *_taboo ) { taboos = _taboos; }
            //! Sets objective pointer
            set(  typename t_ObjectiveType::Vector*  _o )
              { objective = _o; }
            //! Sets objective pointer
            set(  typename t_Store::Base*  _s ) { store = _s; }
            //! Sets history pointer
            set(  typename t_Store::Base*  _history ) { history = _history; }
        };
        
        //! \brief Communication CRT class for bulls.
        //! \details This class defines the requests that a bull can send a
        //!          farmer, and the commands it can send to the cows (as a
        //!          whole). The latter are:
        //!          - Bull::t_CowCommands::CONTINUE, 
        //!          - Bull::t_CowCommands::DONE, 
        //!          - Bull::t_CowCommands::EVALUATE, 
        //!          - Bull::t_CowCommands::GRADIENT, 
        //!          - Bull::t_CowCommands::WITH_GRADIENT, 
        //!          - Bull::t_CowCommands::ONE_GRADIENT.
        //!          .
        //!          A number of helper functions are also declared for
        //!          broadcasting stuff to cows and requesting stuff from the
        //!          farmer.
        template< T_DERIVED >
        class Bull : private ::mpi::Base
        {
          friend class Cow<T_GATRAITS, T_DERIVED>;
          public:
            //! Type of the derived class
            typedef T_DERIVED t_Derived;
        
          protected:
            //! All %GA traits
            typedef typename t_Derived::t_GATraits t_GATraits;
            //! Type of the farmer communication class
            typedef Farmer<t_GATraits, t_Derived> t_Farmer;
            //! Type of individual in this %GA
            typedef typename t_GATraits :: t_Individual         t_Individual;
            //! Type of the fitness, as declared in the base class
            typedef typename t_Base :: t_Fitness                t_Fitness;
            //! Type of the quantity, as declared in the base class
            typedef typename t_Base :: t_Quantity               t_Quantity;
            //! Type of the scalar quantity, as declared in the base class
            typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
            //! Type of the lamarckian traits, as declared in the base class
            typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
            //! Type of the lamarckian gradients, as declared in the base class
            typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
            //! Type of the lamarckian variables, as declared in the base class
            typedef typename t_VA_Traits :: t_Type              t_VA_Type;
            //! Requests to send to farmer.
            typedef t_Farmer::t_Requests t_Requests;
            //! Commands received from farmer.
            typedef t_Farmer::t_Commands t_Commands;
            //! Commands to send to cows.
            enum t_CowCommands 
            {
              CONTINUE, //!< keep doing what you're doing.
              DONE,     //!< stop doing what you're doing.
              EVALUATE, //!< Evaluate and individual.
              WITH_GRADIENT, //!< Evaluate and compute gradient of an individual.
              GRADIENT, //!< Compute gradient of an individual.
              ONE_GRADIENT, //!< Compute gradient in one direction of an individual.
            }
            //! Type of the base class
            typedef ::mpi::Base t_Base;
        
            //! Tag for communication with cows.
            const static MPI::INT COWTAG = 2;
        
          protected:
            MPI::Comm *cowcomm;
        
          protected:
            Bull   ( MPI::Comm *_fcomm, MPI::Comm *_ccom )
                 : t_Base( _fcom ), cowcomm( _ccom ) {}
        
            //! Sends an individual to farmer.
            void send_individual( t_Individual &_indiv)
              { send_template_object< t_Individual >( 0, _indiv, t_CommBase::comm ); }
            //! Receives an individual from farmer.
            void receive_individual( t_Individual &_indiv )
              { receive_template_object< t_Individual >( 0, _indiv, t_CommBase::comm ); }
            //! Sends an fitness to farmer.
            void send_fitness( types::t_unsigned _bull, t_Fitness &_fit )
              { send_template_object< t_Fitness >( 0, _fit, t_CommBase::comm ); }
            //! Receives a fitness from farmer.
            void receive_fitness( t_Fitness &_fit )
              { receive_template_object< t_Fitness >( 0, _fit, t_CommBase::comm ); }
            //! Sends quantities to farmer.
            void send_quantities( t_Quantities &_q );
              { send_quantities< t_Quantities >( 0, _q, t_CommBase::comm ); }
            //! Sends gradients to farmer.
            void send_gradients( t_QuantityGradients &_grad );
              { send_quantities< t_QuantityGradients >( 0, _grad, t_CommBase::comm ); }
            //! Receives gradients from farmer
            void receive_gradients( t_VA_type *_grad, types::t_unsigned _n );
              { receive_range< t_VA_type >( 0, _grad, *_grad + _n, t_CommBase::comm ); }
            //! Sends an individual to bull \a _bull.
            void receive_quantities( t_Quantities &_fit );
              { receive_quantities< t_Fitness >( 0, _q, t_CommBase::comm ); }
            //! Sends an object to \a _bull.
            template < class t_OBJECT > void send_object( T_OBJECT _object );
              { send_object< T_OBJECT >( 0, _object, t_CommBase::comm ); }
            //! Receives an object from \a _bull.
            template < class t_OBJECT > void receive_object( T_OBJECT &_object );
              { receive_object< T_OBJECT >( 0, _object, t_CommBase::comm ); }
            //! Receives a command from Farmer.
            t_Commands obey();
            //! Sends a request to Farmer.
            void request( t_Requests _request );
            //! Broadcasts a command to all cows
            void command( t_Commands _c );
            //! Broadcasts an individual to all cows
            void bcast( t_Individual &_individual )
              { bcast_template_object( 0, _indiv, _comm ); }
        };
        
        
        //! \brief Communication CRT class for cows.
        //! \details Defines the behavior of a cow responding to commands from
        //!          its bull. Although is meant to act a CRT, and since the
        //!          behavior of cows is expected not to vary much, the
        //!          complete behaviors are defined here. As such, it is
        //!          important to call the routine Cow::set() prior to use.
        template< T_DERIVED >
        class Cow : private ::mpi::Base
        {
          public:
            //! Type of the derived class
            typedef T_DERIVED t_Derived;
        
          private:
            //! All %GA traits
            typedef typename t_Derived::t_GATraits             t_GATraits;
            //! Type of the evaluator
            typedef typename T_GATRAITS :: t_Evaluator         t_Evaluator;
            //! \brief quantity traits pertaining to Virtual Atom minimization
            typedef typename t_GATraits :: t_VA_Traits         t_VATraits;
            //! \brief Gradients type for Virtual Atom minimization
            typedef typename t_VATraits :: t_QuantityGradients t_QuantityGradients;
        
          protected:
            //! Type of the bull communication class
            typedef Bull<t_GATraits, t_Derived> t_Bull;
            //! Commands received from farmer.
            typedef t_Bull::t_Cowands t_Commands;
            //! Type of the base class
            typedef ::mpi::Base t_Base;
            //! Gradients for minimization
            t_QuantityGradients gradients;
            //! Tag for communication with bull.
            const static MPI::INT TAG = t_Bull :: COWTAG;
        
          protected:
            //! Pointer to the interface to the functional(s).
            t_Evaluator *evaluator;
        
          protected:
            //! Constructor and Initializer
            Cow( MPI::Comm *_comm ): t_Base( _comm ), evaluator(NULL) {}
        
            //! Wait for a command from the bull
            t_Commands obey();
            //! Evaluates an individual.
            void onEvaluate();
            //! Evaluates the gradient of an individual.
            void onGradient();
            //! Evaluates an individual and its gradient.
            void onWithGradient();
            //! Evaluates the gradient of an individual in one direction.
            void onOneGradient();

          public:
            //! Sets the pointer to the evaluator
            void set( t_Evaluator *_eval ) { evaluator = _eval; }
        };

        //! \brief Your average cow class.
        //! \details In practive, a cow waits for orders while it watches the
        //!          trains go by. As such its base behavior can be defined quite easily.
        template< class T_GATRAITS, template < class > class T_BASE>
        class LaNormande : private Comm::Bull< Farmer<T_BASE> >,
                           public T_BASE< T_GATRAITS >
        {
          public:
            //! Base class type with meta-evaluator
            typedef T_BASE<t_GATraits> t_Base;
            //! all %GA traits
            typedef typename t_Base::t_GATraits t_GATraits;
    
          protected:
            //! all individual traits
            typedef typename t_GATraits::t_IndivTraits t_IndivTraits;
            //! type of an individual
            typedef typename t_GATraits::t_Individual  t_Individual; 
            //! type of the population
            typedef typename t_GATraits::t_Population  t_Population; 
            //! Communication base class
            typedef Comm::Cow< T_GATRAITS, Cow > t_CommBase;
    
          public:
            //! Constructor.
            LaNormande   ( Topology *_topo )
                       : t_CommBase(_topo ), t_Base() {}
    
            //! Creates \a _offspring population from \a _parent
            void operator()(const t_Population& _parents, t_Population& _offspring);
              { while( t_CommBase :: wait() != t_CommBase::DONE ); }
       
            //! The class name. EO required
            virtual std::string className() const
              { return "GA::mpi::Graph::Comm::LaNormande"; }
        };
      } // namespace Comm
    } // namespace Graph

  } // namespace mpi
} // namespace GA

#include "comm.impl.h"

#endif
#endif
