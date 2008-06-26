//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctime>

#include "random.h"

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>


#include "debug.h"


namespace opt
{
  namespace random
  {
    namespace details
    {

      //! A random number generator.
      boost::mt11213b *generator = NULL;
      //! A uniform distribution.
      boost::uniform_real<> *uni_dist = NULL;
      //! A random number generator
      boost::variate_generator< boost::mt11213b&, 
                                        boost::uniform_real<> > *rng = NULL;
    }

    void create( types::t_unsigned _seed )
    {
      details::generator = new boost::mt11213b(42u);
      if( not _seed ) details::generator->seed( static_cast<unsigned int>( std::time(0) ) );
      else details::generator->seed( static_cast<unsigned int>(_seed) );

      details::uni_dist = new boost::uniform_real<>(0,1);

      details::rng = new boost::variate_generator< boost::mt11213b&, 
                                                   boost::uniform_real<> > 
                                                 ( *details::generator, *details::uni_dist );

      seed( _seed );
    }

    void seed( types::t_unsigned _seed )
    {
      __ASSERT( not details::generator, "random number generator was not created" );
      if( not _seed ) details::generator->seed( static_cast<unsigned int>( std::time(0) ) );
      else details::generator->seed( static_cast<unsigned int>(_seed) );
    }

    types::t_real rng() 
    {
      __ASSERT( not details::rng, "random number generator was not created" );
      return (*details::rng)();
    } 

    void destroy()
    {
      if( details::uni_dist ) delete details::uni_dist;
      details::uni_dist = NULL;
      if( details::generator ) delete details::generator;
      details::generator = NULL;
      if( details::rng ) delete details::rng;
      details::rng = NULL;
    }
  }
}
