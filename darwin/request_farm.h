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
    class Farmer : private ::mpi::Base
    {
      protected:
        struct t_ToBull
        {
          const static ::mpi::INT UNDEFINED = -1;
          const static ::mpi::INT EXITLOOP = 0;
          const static ::mpi::INT EVALUATEINDIVIDUAL = 1;
        };
        struct t_FromBull
        {
          const static ::mpi::INT UNDEFINED = -1;
          const static ::mpi::INT WAITING = 0;
          const static ::mpi::INT SENDINGRESULT = 1;
          const static ::mpi::INT REQUESTINGTABOOCHECK = 2;
        };
        const static types::t_int TAG = 1;
      
      protected:
        MPI::Prequest *to_bulls;
        MPI::Prequest *from_bulls;
        types::t_unsigned *in;
        types::t_unsigned *out;
        types::t_unsigned nbulls;


      public:
        Farmer( MPI :: Comm &_comm );
        ~Farmer();
    
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
