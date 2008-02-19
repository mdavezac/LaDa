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
    template< T_GATRAITS >
    class FarmerComm : private ::mpi::Base
    {
      public:
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
      protected:
        struct t_Operations
        {
          const static ::mpi::INT WAITING = 0;
          const static ::mpi::INT REQUESTINGOBJECTIVE = 1;
          const static ::mpi::INT REQUESTINGTABOOCHECK = 2;
          const static ::mpi::INT REQUESTINGHISTORYCHECK = 2;
        };
        const static types::t_int TAG = 1;
      
      protected:
        //! Holds requests from bull;
        MPI::Prequest *from_bulls;
        //! Request buffers
        types::t_unsigned *requests;
        //! Taboo functor.
        Taboo_Base<t_Individual>*          taboos;
        //! Objective functor.
        typename t_ObjectiveType::Vector*  objective;
        //! Store functor.
        typename t_Store :: Base*          store;


      public:
        Farmer( MPI :: Comm &_comm );
        ~Farmer();

        bool ProbeBulls()
        types::t_real ReceiveReal( types::t_unsigned _bull );
    
    };

    class Bull : private ::mpi::Base
    {

      protected:
        struct t_FromFarmer
        {
          const static ::mpi::INT EXITLOOP = 0;
          const static ::mpi::INT EVALUATEINDIVIDUAL = 1;
        };
        struct t_ToFarmer
        {
          const static ::mpi::INT WAITING = 0;
          const static ::mpi::INT SENDINGRESULT = 1;
          const static ::mpi::INT REQUESTINGTABOOCHECK = 2;
        };
        struct t_ToCows
        {
          const static ::mpi::INT KEEPGOING = 0;
          const static ::mpi::INT EXITLOOP = 1;
          const static ::mpi::INT RECEIVEINDIVIDUAL = 2;
        }

      protected:
        MPI::Comm *cowcomm;
        MPI::Request to_farmer;
        MPI::Request from_farmer;
        types::t_unsigned farmer_in;
        types::t_unsigned farmer_out;

        const static types::t_int COWTAG = 1;
        const static types::t_int BULLTAG = 2;
        
      public:
        Bull();
    
    };

    class Cow : private ::mpi::Base
    {
    };
  }
} 

#endif
#endif
