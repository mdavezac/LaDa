//
//  Version: $Id$
//
#ifdef _MPI
#ifndef _MPI_REQUESTFARM_H_
#define _MPI_REQUESTFARM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <opt/types.h>
#include <mpi/mpi_object.h>

namespace BitString
{
  //! MPI for bitstrings
  namespace mpi
  {
    //! Header to a general bitstring of integers, reals, and characters
    struct Header
    {
      //! size of the array of reals
      types::t_unsigned realsize;
      //! size of the array of integers
      types::t_unsigned intsize;
      //! size of the array of characters
      types::t_unsigned charsize;
    };
    //! Bitstring of integers, reals, and characters
    struct MPIData
    {
      //! array of reals
      types::t_real *realbits;
      //! array of integers
      types::t_int *intbits;
      //! array of characters
      types::t_char *charbits;
    }
  }
}

namespace GA
{

  //! MPI for the Genetic Algorithm.
  namespace mpi
  {
    //! Sends a template object to bull \a _bull.
    template< class T_QUANTITY >
    void send_quantity( types::t_unsigned _bull, T_QUANTITY &_q, MPI::Comm *_comm);
    //! Receives a template object from bull \a _bull.
    template< class T_QUANTITY >
    void receive_quantity( types::t_unsigned _bull, T_QUANTITY &_q, MPI::Comm *_comm);
    //! Sends a template object to bull \a _bull.
    template< class T_OBJECT >
    void send_template_object( types::t_unsigned _bull, T_OBJECT &_o, MPI::Comm *_comm);
    //! Receives a template object from bull \a _bull.
    template< class T_OBJECT >
    void receive_template_object( types::t_unsigned _bull, T_OBJECT &_o, MPI::Comm *_comm);
    //! Sends a (non-template) object to bull \a _bull.
    template< class T_OBJECT >
    void send_template_object( types::t_unsigned _bull, T_OBJECT &_o, MPI::Comm *_comm);
    //! Receives a (non-template) object from bull \a _bull.
    template< class T_OBJECT >
    void receive_template_object( types::t_unsigned _bull, T_OBJECT &_o, MPI::Comm *_comm);
    //! Broadcasts a template object from \a _root.
    template< class T_OBJECT >  
    void bcast_template_object( types::t_int _root, T_OBJECT &_object, MPI::Comm *_comm )

    template< T_GATRAITS, T_DERIVED > class CowComm;
    template< T_GATRAITS, T_DERIVED > class BullComm;

    template< T_GATRAITS, T_DERIVED >
    class FarmerComm : private ::mpi::Base
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
          REQUESTINGOBJECTIVE, //!< Requesting an objective evaluation.
          REQUESTINGTABOOCHECK, //!< Requesting to know whether an individual is taboo.
          REQUESTINGHISTORYCHECK; //!< Requesting to know whether an individual is known.
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

      public:
        FarmerComm( MPI :: Comm *_comm );
        ~FarmerComm();

        //! \brief Tests incoming messages from bull and dispatches response.
        //! \details The incoming messages are listed in
        //!          FarmerComm::t_Requests. There response are dispatched to
        //!          to onWait(), onTaboo(), onHistory(), and onCheck() through
        //!          a CRT mechanism (Curriuosly Recurring Template), e.g. to
        //!          the member functions of the derivatives of this class.
        bool test_bulls();

        //! Sets taboo pointer
        set( Taboo_Base<t_Individual> *_taboo ) { taboos = _taboos; }
        //! Sets objective pointer
        set(  typename t_ObjectiveType::Vector*  _o )
          { objective = _o; }
        //! Sets objective pointer
        set(  typename t_Store::Base*  _s ) { store = _s; }
        //! Sets history pointer
        set(  typename t_Store::Base*  _history ) { history = _history; }

        //! Sends an individual to bull \a _bull.
        void send_individual( types::t_unsigned _bull, t_Individual &_indiv)
          { send_template_object< t_Individual >( _bull, _indiv, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void receive_individual( types::t_unsigned _bull, t_Individual &_indiv )
          { receive_template_object< t_Individual >( _bull, _indiv, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void send_fitness( types::t_unsigned _bull, t_Fitness &_fit )
          { send_template_object< t_Fitness >( _bull, _fit, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void receive_individual( types::t_unsigned _bull, t_Fitness &_fit )
          { receive_template_object< t_Fitness >( _bull, _fit, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void send_quantities( types::t_unsigned _bull, t_Quantities &_q );
          { send_quantities< t_Fitness >( _bull, _q, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void receive_quantities( types::t_unsigned _bull, t_Quantities &_fit );
          { receive_quantities< t_Fitness >( _bull, _q, t_CommBase::comm ); }
        //! Starts all persistent requests from bulls ( FarmerComm::requests )
        void startall() { _comm->StartAll( nbulls, requests ); } 
        //! Sends a command to \a _bull.
        void send_command( types::t_unsigned _bull, t_Commands _c );
    };

    template< T_GATRAITS, T_DERIVED >
    class BullComm : private ::mpi::Base
    {
      friend class CowComm<T_GATRAITS, T_DERIVED>;
      public:
        //! Type of the derived class
        typedef T_DERIVED t_Derived;
        //! All %GA traits
        typedef typename T_GATRAITS t_GATraits;

      protected:
        //! Type of the farmer communication class
        typedef FarmerComm<t_GATraits, t_Derived> t_FarmerComm;
        //! Requests to send to farmer.
        typedef t_FarmerComm::t_Requests t_Requests;
        //! Commands received from farmer.
        typedef t_FarmerComm::t_Commands t_Commands;
        //! Commands to send to cows.
        enum t_CowCommands 
        {
          CONTINUE, //!< keep doing what you're doing.
          DONE,     //!< stop doing what you're doing.
          EVALUATE, //!< Evaluate and individual.
          EVALUATE_WITH_GRADIENT, //!< Evaluate and compute gradient of an individual.
          EVALUATE_GRADIENT, //!< Compute gradient of an individual.
          EVALUATE_ONE_GRADIENT, //!< Compute gradient in one direction of an individual.
        }
        //! Type of the base class
        typedef ::mpi::Base t_Base;

        //! Tag for communication with cows.
        const static MPI::INT COWTAG = 2;

      protected:
        MPI::Comm *cowcomm;

      public:
        BullComm( MPI::Comm *_fcomm, MPI::Comm *_ccom ) : t_Base( _fcom ), cowcomm( _ccom ) {}

        //! Sends an individual to bull \a _bull.
        void send_individual( t_Individual &_indiv)
          { send_template_object< t_Individual >( 0, _indiv, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void receive_individual( t_Individual &_indiv )
          { receive_template_object< t_Individual >( 0, _indiv, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void send_fitness( types::t_unsigned _bull, t_Fitness &_fit )
          { send_template_object< t_Fitness >( 0, _fit, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void receive_individual( t_Fitness &_fit )
          { receive_template_object< t_Fitness >( 0, _fit, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void send_quantities( t_Quantities &_q );
          { send_quantities< t_Fitness >( 0, _q, t_CommBase::comm ); }
        //! Sends an individual to bull \a _bull.
        void receive_quantities( t_Quantities &_fit );
          { receive_quantities< t_Fitness >( 0, _q, t_CommBase::comm ); }
        //! Receives a command from Farmer.
        t_Commands obey();
        //! Sends a request to Farmer.
        void request( t_Requests _request );
        //! Broadcasts a command to all cows
        void bcast( t_Commands _c );
        //! Broadcasts and individual to all cows
        void bcast( t_Individual &_individual )
          { bcast_template_object( 0, _indiv, _comm ); }
    
    };


    template< T_GATRAITS, T_DERIVED >
    class CowComm : private ::mpi::Base
    {
      public:
        //! Type of the derived class
        typedef T_DERIVED t_Derived;
        //! All %GA traits
        typedef T_GATRAITS t_GATraits;

      private:
        //! Type of the evaluator
        typedef typename T_GATRAITS :: t_Evaluator t_Evaluator;
        //! \brief quantity traits pertaining to Virtual Atom minimization
        typedef typename t_GATraits :: t_VA_Traits                   t_VATraits;
        //! \brief Gradients type for Virtual Atom minimization
        typedef typename t_VATraits :: t_QuantityGradients           t_QuantityGradients;

      protected:
        //! Type of the bull communication class
        typedef BullComm<t_GATraits, t_Derived> t_BullComm;
        //! Commands received from farmer.
        typedef t_BullComm::t_CowCommands t_Commands;
        //! Type of the base class
        typedef ::mpi::Base t_Base;
        //! Gradients for minimization
        t_QuantityGradients gradients;
        //! Tag for communication with bull.
        const static MPI::INT TAG = t_BullComm :: COWTAG;

      protected:
        t_Evaluator *evaluator;

      public:
        //! Constructor and Initializer
        CowComm( MPI::Comm *_comm ): t_Base( _comm ), evaluator(NULL) {}

        //! Sets the pointer to the evaluator
        void set( t_Evaluator *_eval ) { evaluator = _eval; }

        //! Wait for a command from the bull
        t_Commands obey();
        //! Evaluates an individual.
        void onEvaluate();
    };

  }
} 

#endif
#endif
