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
#include <boost/shared_ptr.hpp>


#include "debug.h"


namespace opt
{
  namespace random
  {
    namespace details
    {

      //! A random number generator.
      boost::shared_ptr< boost::mt11213b > generator;
      //! A uniform distribution.
      boost::shared_ptr< boost::uniform_real<> > uni_dist;
      //! A random number generator
      boost::shared_ptr< boost::variate_generator< boost::mt11213b&, 
                                                   boost::uniform_real<> > > rng;
    }

    void create( types::t_unsigned _seed )
    {
      boost::shared_ptr< boost::mt11213b > a( new boost::mt11213b(42u) );
      __DOASSERT( not a.get(), "Could not create random number generator boost::mt11213b.\n" )
      details::generator = a;
      if( not _seed ) details::generator->seed( static_cast<unsigned int>( std::time(0) ) );
      else details::generator->seed( static_cast<unsigned int>(_seed) );

      boost::shared_ptr< boost::uniform_real<> > b( new boost::uniform_real<>(0,1) );
      __DOASSERT( not a.get(), "Could not create uniform distribution.\n" )
      details::uni_dist = b;
      
      typedef boost::variate_generator< boost::mt11213b&, 
                                        boost::uniform_real<> > t_Shit;
      boost::shared_ptr< t_Shit > c( new t_Shit( *details::generator, *details::uni_dist ) );
      __DOASSERT( not a.get(), "Could not create random number "
                               "generator with uniform distribution.\n" )
      details::rng = c; 
                                                 
      seed( _seed );
    }

    void seed( types::t_unsigned _seed )
    {
      if( not details::generator.get() )
      {
        create( _seed );
        return;
      } 
      __ASSERT( not details::generator.get(), "random number generator was not created" );
      if( not _seed ) details::generator->seed( static_cast<unsigned int>( std::time(0) ) );
      else details::generator->seed( static_cast<unsigned int>(_seed) );
    }

    types::t_real rng() 
    {
      __ASSERT( not details::rng.get(), "random number generator was not created" )
      return (*details::rng)();
    } 
    bool flip()
    {
      __ASSERT( not details::rng.get(), "random number generator was not created" )
      return (*details::rng)() -0.5e0 > 0e0;
    }
    types::t_real rflip()
    {
      __ASSERT( not details::rng.get(), "random number generator was not created" )
      return 2e0*(*details::rng)() - 1e0;
    }

    types::t_unsigned range( types::t_unsigned _first, types::t_unsigned _last )
    {
      __ASSERT( not details::rng.get(), "random number generator was not created" )
      return   (types::t_unsigned ) 
                    std::floor( (*details::rng)() * (types::t_real)( _last - _first ) )
             + _first;
    }
  }
}
